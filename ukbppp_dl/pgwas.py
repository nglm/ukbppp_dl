"""Functions for processing pGWAS summary statistics files from UKB PPP."""

import sys
import time
import synapseclient
from synapseclient.models import Folder
import tarfile
import gzip
import polars as pl
import os
import json
import pathlib
from pathlib import Path
from datetime import datetime
from typing import IO, Any, Dict, List, Optional, Tuple, Union

from .common import save_log, delete_files, download_from_synapse

_CONSISTENCY_MESSAGE = """
    Files were skipped because their result file already exists.
    Please make sure that the pre-existing result files are correct and consistent with the new analysis.
    If not consistent, delete result files and re-run the analysis.
    Decisive parameters include: log10p_threshold, regenie_columns, csv_columns.
"""

PGWAS_REGIONS = {
    "African" : "syn51365304",
    "American" : "syn51500434",
    "Central_South_Asian" : "syn51365305",
    "Combined" : "syn51365308",
    "East_Asian" : "syn51365306",
    "European" : "syn51365303",
    "Middle_East" : "syn51365307",
}

def list_tar_files_in_region_folder(
        synapse_id: str = PGWAS_REGIONS["Combined"],
        login_kwargs: Dict[str, Any] = {},
    ) -> List[Tuple[str, str]]:
    """
    List all protein tar files in the folder of given GWAS region.

    Parameters
    ----------
    synapse_id : str, optional
        Synapse folder ID for the GWAS region, e.g. ``"syn51365308"``
        for the Combined ancestry group. Defaults to
        ``PGWAS_REGIONS["Combined"]``.
    login_kwargs : Dict[str, Any], optional
        Keyword arguments forwarded to
        ``synapseclient.Synapse.login()``. Common keys include
        ``"authToken"`` (``str``), ``"email"`` (``str``) and
        ``"profile"`` (``str``). See `synapseclient.Synapse.login()
        <https://python-docs.synapse.org/en/stable/reference/client/#synapseclient.Synapse.login>`_
        for more information.

        Alternatively, you can configure your Synapse credentials in a
        ``.synapseConfig`` file and leave this argument empty to use the
        default login behaviour. You can find more information about
        ``.synapseConfig`` in the `Synapse documentation
        <https://python-docs.synapse.org/en/stable/tutorials/authentication/>`_.

    Returns
    -------
    List[Tuple[str, str]]
        Sorted list of ``(synapse_id, tar_filename)`` pairs, one per
        protein tar file found in the folder. Example entry:
        ``[("syn52361344",
        "ABCA2_Q9BZC7_OID30146_v1_Cardiometabolic_II.tar"),
        ("syn51470065",
        "ABHD14B_Q96IU4_OID20921_v1_Neurology.tar.tar")]``.
    """
    syn = synapseclient.Synapse()
    # Use your synapse account here, as configured in your .synapseConfig file.
    syn.login(**login_kwargs)

    folder = Folder(id=synapse_id)
    # Iterator of (dirpath, dirs, nondirs)
    folder_structure = folder.walk(recursive=False)

    # There is actually only one level of files in the folder,
    # So we directly get the first iteration of the walk
    # nondirs is a list of EntityHeader objects, which have attributes such as id and name
    dirpath, dirs, nondirs = next(folder_structure)

    # Filter for tar files and return their Synapse IDs and names
    tar_files = [
        (file.id, file.name)
        for file in nondirs if file.name.endswith('.tar')
    ]

    tar_files.sort(key=lambda x: x[1])  # Sort by file name

    return tar_files


def keep_significant_qtls_from_chr_gz_file(
        chr_file_gz: str | os.PathLike[str] | IO[bytes],
        separator: str = " ",
        columns: list[str] | None = None,
        new_columns: list[str] | None = None,
        log10p_threshold: float | int = 7,
        create_log: bool | int = False,
        log_kwargs: dict[str, Any] = {},
        verbose: bool | int = False,
    ) -> tuple[pl.DataFrame, dict]:
    """
    Keep significant QTLs from a gzipped summary-statistics file.

    May create a JSON log file.

    Parameters
    ----------
    chr_file_gz : str or os.PathLike[str] or IO[bytes]
        Path to the ``.gz`` file, or a binary file-like object (e.g.
        extracted from a ``tarfile``). To be passed directly to
        ``gzip.open``.
    separator : str, optional
        Field separator used in the REGENIE file. Default is ``" "``
        (space).
    columns : list[str] or None, optional
        Subset of column names to load from the REGENIE file, e.g.
        ``["CHROM", "GENPOS", "ID", "LOG10P"]``. ``None`` loads all
        columns. ``LOG10P`` column must be present.
    new_columns : list[str] or None, optional
        Rename the loaded columns using this list of new names. Must
        have the same length as *columns* when provided. ``LOG10P``
        column must be present and with the same name for filtering to
        work.
    log10p_threshold : float or int, optional
        Minimum ``LOG10P`` value for a QTL to be considered significant.
        Default is ``7`` (roughly p < 1e-7).
    create_log : bool or int, optional
        Whether to save a JSON log file. ``0``/``False`` disables
        logging; ``True`` or any positive integer enables it. Default is
        ``False``.
    log_kwargs : dict[str, Any], optional
        Keyword arguments forwarded to :func:`save_log`. Recognised key:
        ``"log_fname"`` (``str``) to override the default log filename.
    verbose : bool or int, optional
        Verbosity level. ``0``/``False`` is silent; higher values print
        more.

    Returns
    -------
    summary_stats_significant : pl.DataFrame
        DataFrame containing only rows where ``LOG10P >=
        log10p_threshold``.
    log : dict
        Dictionary with keys: ``"log_filename"`` (``str | None``),
        ``"log10p_threshold"`` (``float``), ``"n_tot_qtls"`` (``int``),
        ``"n_kept_qtls"`` (``int``), ``"source_chr_file"`` (``str``).
    """
    # Opening the extracted .gz file to read its content
    with gzip.open(chr_file_gz, 'rt') as f_regenie:

        # Reading the .regenie file as a CSV and creating a dataframe
        summary_stats = pl.read_csv(
            f_regenie,
            separator=separator,
            columns=columns,
            new_columns=new_columns,
            )


        n_tot_qtls = len(summary_stats)

        # Filtering the summary statistics to keep only significant QTLs
        summary_stats_significant = summary_stats.filter(
            pl.col("LOG10P") >= log10p_threshold
        )

        n_kept_qtls = len(summary_stats_significant)

        log = {
            "log_filename": None,
            "log10p_threshold": log10p_threshold,
            "n_tot_qtls": n_tot_qtls,
            "n_kept_qtls": n_kept_qtls,
            "source_chr_file": f_regenie.name,
        }

    if int(create_log) > 0:
        # If log_fname not in log_kwargs, use a default one
        log_fname = f"{Path(chr_file_gz.name).stem}-log.json"
        log_fname = log_kwargs.pop("log_fname", log_fname)
        save_log(log_fname, log, **log_kwargs)

    if int(verbose) > 0:
        print(f"Total QTLs in chr file:    {n_tot_qtls}")
        print(f"Kept QTLs (log10p>= {log10p_threshold}):    {n_kept_qtls}")


    return summary_stats_significant, log

def process_one_chr_from_protein_tar_file(
        protein_tar_file: tarfile.TarFile,
        chr_gz_fname: str,
        res_location: str = "./ukb_ppp_dl/results",
        separator: str = " ",
        columns: list[str] | None = None,
        new_columns: list[str] | None = None,
        log10p_threshold: float | int = 7,
        create_log: bool | int = False,
        log_kwargs: dict[str, Any] = {},
        verbose: bool | int = False,
    ) -> tuple[str, dict]:
    """
    Keep significant QTLs from one chromosome file of a protein tar.

    May create a JSON log file and will create a CSV file.

    Parameters
    ----------
    protein_tar_file : tarfile.TarFile
        Open ``TarFile`` object for the protein's tar archive.
    chr_gz_fname : str
        Name or partial path of the chromosome ``.gz`` file inside the
        tar, e.g. ``"chr1.regenie.gz"``. A substring match is used to
        locate the actual entry.
    res_location : str, optional
        Directory where the result CSV (and log) will be written.
    separator : str, optional
        Field separator used in the REGENIE file. Default is ``" "``
        (space).
    columns : list[str] or None, optional
        Subset of column names to load from the REGENIE file, e.g.
        ``["CHROM", "GENPOS", "ID", "LOG10P"]``. ``None`` loads all
        columns. ``LOG10P`` column must be present.
    new_columns : list[str] or None, optional
        Rename the loaded columns using this list of new names. Must
        have the same length as *columns* when provided. ``LOG10P``
        column must be present and with the same name for filtering to
        work.
    log10p_threshold : float or int, optional
        Minimum ``LOG10P`` threshold for significance. Default is ``7``.
    create_log : bool or int, optional
        Verbosity of per-chromosome logging. ``0`` disables; ``1``
        creates a JSON log for this chromosome file.
    log_kwargs : dict[str, Any], optional
        Keyword arguments forwarded to :func:`save_log`. Recognised key:
        ``"log_fname"`` (``str``) to override the default log filename.
    verbose : bool or int, optional
        Verbosity level.

    Returns
    -------
    res_csv_fname : str
        Path to the CSV file containing significant QTLs for this
        chromosome.
    log_chr : dict
        Log dictionary with key ``"skipped"`` (``bool``), and when not
        skipped: ``"log10p_threshold"`` (``float``), ``"n_tot_qtls"``
        (``int``), ``"n_kept_qtls"`` (``int``), ``"source_chr_file"``
        (``str``).
    """

    # Find actual full name in case chr_gz_fname doesn't contain the full path in the tar file
    list_of_chr_files = protein_tar_file.getnames()
    chr_gz_fname_full = [f for f in list_of_chr_files if chr_gz_fname in f][0]
    chr_gz_fname_no_path_no_gz = f"{Path(chr_gz_fname_full).stem}"

    # Preparing result filename keeping only significant QTLs
    res_csv_fname = f"{res_location}/{chr_gz_fname_no_path_no_gz}_significant_qtls.csv"
    p = pathlib.Path(res_csv_fname)
    p.parent.mkdir(parents=True, exist_ok=True)

    # Don't run again if result file already exists
    if os.path.isfile(f"{res_csv_fname}"):

        log_chr = {"skipped": True}

        if int(verbose) > 0:
            print(
                f"\n\tResults for chromosome file {chr_gz_fname} already created at {res_csv_fname}. Skipping.",
                flush=True
            )

    else:

        if int(verbose) > 0:
            print(f"\n\tProcessing chromosome file:     {chr_gz_fname}.")

        t_start = time.time()

        # Temporarily extracting the .gz file to read its content
        chr_file_gz = protein_tar_file.extractfile(chr_gz_fname_full)

        # Read .regenie_file and keep only significant QTLs
        summary_stats_significant, log_chr = keep_significant_qtls_from_chr_gz_file(
            chr_file_gz=chr_file_gz,
            separator=separator,
            columns=columns,
            new_columns=new_columns,
            log10p_threshold=log10p_threshold,
            create_log=False,  # create the log after merge
            log_kwargs=log_kwargs,
            verbose=int(verbose) - 1,
        )

        log_chr["skipped"] = False

        # Saving the significant summary statistics to a new file
        summary_stats_significant.write_csv(f"{res_csv_fname}")
        t_end = time.time()

        if int(create_log) > 0:
            # If log_fname not in log_kwargs, use a default one
            log_fname = f"{res_location}/{chr_gz_fname_no_path_no_gz}-log.json"
            log_fname = log_kwargs.pop("log_fname", log_fname)
            save_log(log_fname, log_chr, **log_kwargs)

        if int(verbose) > 0:
            print(f"\tSignificant QTLs saved to:      {res_csv_fname}", flush=True)
            print(f"\tProcessed chromosome file in:   {t_end - t_start:.2f} s", flush=True)

    return res_csv_fname, log_chr


def merge_significant_qtls_from_csv(
        csv_fnames: list[str],
        output_fname: str | None = None,
        create_log: bool | int = False,
        log_kwargs: dict[str, Any] = {},
        delete_csv: bool = False,
        verbose: bool | int = False,
    ) -> tuple[pl.DataFrame, dict]:
    """
    Concatenate multiple CSV files into one DataFrame.

    May create a JSON log file, may create a CSV file and can optionally
    delete the input CSV files after merging.

    Parameters
    ----------
    csv_fnames : list[str]
        Paths to the per-chromosome CSV files to merge. Must be
        non-empty. Each file must contain at least a ``LOG10P`` column.
        Files must be consistent with each other (e.g. same columns) but
        can have no rows (but a header).
    output_fname : str or None, optional
        If provided, the merged DataFrame is written to this CSV path.
    create_log : bool or int, optional
        If truthy *and* *output_fname* is set, save a JSON log alongside
        the merged CSV.
    log_kwargs : dict[str, Any], optional
        Keyword arguments forwarded to :func:`save_log`. Recognised key:
        ``"log_fname"`` (``str``) to override the default log filename.
    delete_csv : bool, optional
        If ``True``, delete the input *csv_fnames* after merging.
        Default is ``False``.
    verbose : bool or int, optional
        Verbosity level.

    Returns
    -------
    all_significant_qtls : pl.DataFrame
        Concatenated DataFrame of all significant QTLs from every input
        file. Empty (but with a header) if no significant QTL was
        found.
    log : dict
        Dictionary with keys: ``"log_filename"`` (``str | None``),
        ``"merged_csv_filename"`` (``str | None``), ``"n_csv_merged"``
        (``int``), ``"n_kept_qtls"`` (``int``),
        ``"n_kept_qtls_per_csv"`` (``dict[str, int]``, mapping each CSV
        path to its row count), ``"min_log10p"`` (``float | None``).
    """
    # Note: LOG10P column must exist in the csv files and with the correct case

    # To avoid a too different case that should always be handled separately.
    # Having no csv_fnames shouldn't happen anyway because if there is
    # no significant QTL, a empty csv file is still created with a header.
    if not csv_fnames:
        raise ValueError(f"No CSV files given for merging: {csv_fnames}")

    if verbose:
        print(f"    {"-"*50}  ")
        print(f"Merging significant QTLs from CSV files.", flush=True)

    # Load significant QTLs from the csv files and concatenate them
    all_df = [pl.read_csv(csv_fname) for csv_fname in csv_fnames]
    n_kept_per_file = [len(df) for df in all_df]

    # Keep only non-empty dataframes to concatenate and avoid error
    non_empty_df = [df for df in all_df if len(df) > 0]

    # If there is no non-empty dataframe, create an empty dataframe with the same columns as the input dataframes, to avoid errors later on
    if len(non_empty_df) > 0:
        all_significant_qtls = pl.concat(non_empty_df)
    else:
        # Create an empty df with the same columns as the input dataframes
        all_significant_qtls = all_df[0].clear()

    n_kept_qtls = len(all_significant_qtls)
    if n_kept_qtls > 0:
        min_log10p = float(all_significant_qtls.get_column("LOG10P").min())
    else:
        min_log10p = None

    log = {
        "log_filename" : None,
        "merged_csv_filename" : None,
        "n_csv_merged": len(csv_fnames),
        "n_kept_qtls": n_kept_qtls,
        "n_kept_qtls_per_csv": dict(zip(csv_fnames, n_kept_per_file)),
        "min_log10p": min_log10p,
    }

    # Write the merged significant QTLs to csv if output_fname is specified
    if output_fname is not None:
        # Save the concatenated dataframe to a new file
        all_significant_qtls.write_csv(output_fname)
        log["merged_csv_filename"] = output_fname

    # Save log only if we also create a output file
    if int(create_log) > 0 and output_fname is not None:
        # If log_fname not in log_kwargs, use a default one
        log_fname = f"{output_fname.split('.csv')[0]}-log.json"
        log_fname = log_kwargs.pop("log_fname", log_fname)
        save_log(log_fname, log, **log_kwargs)

    if delete_csv:
        n_deleted = delete_files(csv_fnames, verbose=int(verbose) - 1)
        if int(verbose) > 0 and n_deleted > 0:
            print(f"Deleted {n_deleted}/{len(csv_fnames)} CSV files.")

    if int(verbose) > 0:
        if output_fname:
            print(f"Significant QTLs from {len(csv_fnames)} files saved to:")
            print(f"{output_fname}")
        else:
            print(f"Significant QTLs from {len(csv_fnames)} files concatenated but not saved to file.")
        print(f"Total significant QTLs:         {n_kept_qtls}")
        print(f"  {"-"*50}  \n", flush=True)

    return all_significant_qtls, log

def process_one_tar_file(
        tar_fname: str,
        res_location: str = "./ukb_ppp_dl/results",
        separator: str = " ",
        columns: list[str] | None = None,
        new_columns: list[str] | None = None,
        log10p_threshold: float | int = 7,
        create_log: bool | int = True,
        log_kwargs: dict[str, Any] = {},
        verbose: bool | int = False,
    ) -> tuple[list[str], dict]:
    """
    Keep significant QTLs from one protein tar file.

    May create a JSON log file for the protein and may create JSON log
    files for each chromosome. Will create a CSV file for each
    chromosome.

    Parameters
    ----------
    tar_fname : str
        Path to the protein ``.tar`` file, e.g.
        ``"ABCA2_Q9BZC7_OID30146_v1_Cardiometabolic_II.tar"``.
    res_location : str, optional
        Directory where result CSVs and logs are written.
    separator : str, optional
        Field separator used in REGENIE files. Default is ``" "``
        (space).
    columns : list[str] or None, optional
        Subset of column names to load from the REGENIE file, e.g.
        ``["CHROM", "GENPOS", "ID", "LOG10P"]``. ``None`` loads all
        columns. ``LOG10P`` column must be present.
    new_columns : list[str] or None, optional
        Rename the loaded columns using this list of new names. Must
        have the same length as *columns* when provided. ``LOG10P``
        column must be present and with the same name for filtering to
        work.
    log10p_threshold : float or int, optional
        Minimum ``LOG10P`` threshold for significance. Default is ``7``.
    create_log : bool or int, optional
        Logging verbosity. ``0`` disables; ``1`` creates a tar-level
        JSON log; ``2`` also creates per-chromosome logs. Default is
        ``True`` (``1``).
    log_kwargs : dict[str, Any], optional
        Keyword arguments forwarded to :func:`save_log`. Recognised key:
        ``"log_fname"`` (``str``) to override the default log filename.
    verbose : bool or int, optional
        Verbosity level.

    Returns
    -------
    all_csv_fnames : list[str]
        Paths to the per-chromosome CSV files of significant QTLs.
    log_tar : dict
        Log dictionary with keys: ``"log_filename"`` (``str | None``),
        ``"tar_fname"`` (``str``), ``"protein_name"`` (``str``),
        ``"n_chr_files"`` (``int``), ``"skipped_chr_files"``
        (``list[str]``), ``"all_csv_fnames"`` (``list[str]``),
        ``"n_processed_qtls"`` (``int``), ``"log10p_threshold"``
        (``float``), ``"regenie_columns"`` (``list[str] | None``),
        ``"csv_columns"`` (``list[str] | None``).
    """
    if int(verbose) > 0:
        print(f"Processing protein file: {tar_fname}")

    # There is one .tar file per protein
    with tarfile.open(tar_fname) as protein_tf:

        # In each .tar file, there is one .gz file per chromosome
        list_of_chr_files = protein_tf.getnames()

        # Keep only actual gz files, not the directories
        list_of_chr_files = [f for f in list_of_chr_files if f.endswith('.gz')]
        list_of_chr_files.sort()

        protein_name = f"{Path(tar_fname).stem.split('_')[0]}"

        log_tar = {
            "log_filename": None,
            "tar_fname": tar_fname,
            "protein_name": protein_name,
            "n_chr_files": len(list_of_chr_files),
            "skipped_chr_files": [],
            "all_csv_fnames": [],
            "n_processed_qtls": 0,
            "log10p_threshold": log10p_threshold,
            "regenie_columns": columns,
            "csv_columns": new_columns,
        }

        if verbose:
            print(f"Going through chromosome files...")
        # Focusing on one chromosome at a time
        # For each chr, there one .gz file containing one .regenie file
        for chr_gz_fname in list_of_chr_files:

            res_csv_fname, log_chr = process_one_chr_from_protein_tar_file(
                protein_tar_file=protein_tf,
                chr_gz_fname=chr_gz_fname,
                res_location=res_location,
                separator=separator,
                columns=columns,
                new_columns=new_columns,
                log10p_threshold=log10p_threshold,
                create_log=int(create_log) - 1,
                log_kwargs=log_kwargs,
                verbose=int(verbose) - 1,
            )

            if log_chr["skipped"]:
                log_tar["skipped_chr_files"].append(chr_gz_fname)
            else:
                log_tar["n_processed_qtls"] += log_chr["n_tot_qtls"]
            log_tar["all_csv_fnames"].append(res_csv_fname)


    if log_tar["skipped_chr_files"] and int(verbose) > 0:
        print(f"[WARNING]: Skipped {len(log_tar['skipped_chr_files'])} chromosome files in tar file {tar_fname}!")
        print(_CONSISTENCY_MESSAGE)

    if int(create_log) > 0:
        # If log_fname not in log_kwargs, use a default one
        log_fname = f"{res_location}/{Path(tar_fname).stem}-log.json"
        log_fname = log_kwargs.pop("log_fname", log_fname)
        save_log(log_fname, log_tar, **log_kwargs)

    return log_tar["all_csv_fnames"], log_tar

def find_partial_region_logs(
    example_log_dict: dict | None = None,
    synapse_folder_id: str = PGWAS_REGIONS["Combined"],
    res_location: str = './ukb_ppp_dl/results',
    regenie_columns: list[str] | None = None,
    csv_columns: list[str] | None = None,
    log10p_threshold: float | int = 7,
    all_tar_files: list[str] = [],
    verbose: bool | int = False,
) -> dict[str, dict]:
    """
    Find compatible pre-existing partial region log files.

    Parameters
    ----------
    example_log_dict : dict or None, optional
        If provided, the function reads *synapse_folder_id*,
        *regenie_columns*, *csv_columns*, *log10p_threshold*, and
        *all_tar_files* from this dict instead of the individual keyword
        arguments. Expected keys match those written by
        :func:`keep_significant_qtls_from_region`.
    synapse_folder_id : str, optional
        Synapse folder ID used to filter log filenames, e.g.
        ``"syn51365308"``.
    res_location : str, optional
        Directory to search for ``*.json`` partial log files.
    regenie_columns : list[str] or None, optional
        Expected value of the ``"regenie_columns"`` key in matching logs.
    csv_columns : list[str] or None, optional
        Expected value of the ``"csv_columns"`` key in matching logs.
    log10p_threshold : float or int, optional
        Expected threshold. Logs with a different value are excluded.
    all_tar_files : list[str], optional
        Expected list of tar filenames, e.g.
        ``["ABCA2_Q9BZC7_OID30146_v1_Cardiometabolic_II.tar", ...]``.
    verbose : bool or int, optional
        Verbosity level.

    Returns
    -------
    dict[str, dict]
        Mapping from log file path (``str``) to the parsed log dictionary
        (``dict``) for each compatible partial log found.
    """

    if example_log_dict is not None:
        # If an example log file is provided, we use it to get the parameters to look for in the pre-existing logs

        synapse_folder_id = example_log_dict["synapse_folder_id"]
        regenie_columns = example_log_dict["regenie_columns"]
        csv_columns = example_log_dict["csv_columns"]
        log10p_threshold = example_log_dict["log10p_threshold"]
        all_tar_files = example_log_dict["all_tar_files"]

    # Find all filenames in the result location
    files = os.listdir(res_location)

    # Check that filename corresponds
    log_files = [
        f for f in files
        if f.endswith(".json")
        and f.startswith(f"PART-{synapse_folder_id}-significant_qtl")
    ]
    if verbose:
        print(f"Found {len(log_files)} pre-existing partial log files.")

    compatible_log_files = {}
    non_protein_keys = [
        "synapse_folder_id", "log10p_threshold", "regenie_columns",
        "csv_columns", "all_tar_files", "log_filename", "output_filename",
        "tar_skipped", "tar_processed", "tar_downloaded_skip", "all_proteins",
        "kept_proteins", "partial_log_merged", "partial_output_filenames",
    ]
    to_keep = True

    for log_file in log_files:
        with open(f"{res_location}/{log_file}", 'r') as f_log:
            log_dict = json.load(f_log)

        # ----- Check that essential parameters are consistent -----
        if log_dict["log10p_threshold"] != log10p_threshold:
            to_keep = False
            continue

        if log_dict["regenie_columns"] != regenie_columns:
            to_keep = False
            continue

        if log_dict["csv_columns"] != csv_columns:
            to_keep = False
            continue

        if log_dict["synapse_folder_id"] != synapse_folder_id:
            to_keep = False
            continue

        if set(log_dict["all_tar_files"]) != set(all_tar_files):
            to_keep = False
            continue

        # --------------- Check output text still exist ---------------
        output_fname = log_dict["output_filename"]
        if not os.path.isfile(output_fname):
            to_keep = False
            continue

        # --------------- Check that results still exist ---------------
        protein_keys = [k for k in log_dict.keys() if k not in non_protein_keys]
        for prot in protein_keys:

            res_csv_fname = log_dict[prot]["merged_csv_filename"]
            if not os.path.isfile(res_csv_fname):
                to_keep = False
                continue

        # If we reach this line, it means that the log file is compatible
        if to_keep:
            compatible_log_files[f"{res_location}/{log_file}"] = log_dict

    if int(verbose) > 0:
        print(f"Found {len(compatible_log_files)}/{len(log_files)} partial log files with compatible parameters.")

    return compatible_log_files


def merge_partial_region_logs(
    compatible_log_files: dict[str, dict],
) -> dict:
    """
    Merge multiple compatible partial region log dictionaries into one.

    To be compatible, log files must have been created with the same
    parameters and created by running the function
    :func:`keep_significant_qtls_from_region`. In particular, they must share
    identical ``log10p_threshold``, ``regenie_columns``,
    ``csv_columns``, and ``synapse_folder_id`` values.

    We cannot use this function with already merged log files.

    It is not advised to use this function unless you are extremely sure
    that all log files were created with the same parameters and that
    you will use the merged log file in a consistent way.

    This function is mainly intended to be used as a helper function for
    the :func:`keep_significant_qtls_from_region` function.

    Parameters
    ----------
    compatible_log_files : dict[str, dict]
        Mapping from log file path to log dictionary, as returned by
        :func:`find_partial_region_logs`. All logs must share identical
        ``log10p_threshold``, ``regenie_columns``, ``csv_columns``, and
        ``synapse_folder_id`` values.

    Returns
    -------
    dict:
        A single merged log dictionary. Returns an
        empty ``{}`` when *compatible_log_files* is empty.
    """

    if len(compatible_log_files) == 0:
        return {}
    merged_log = {}
    merging_msg = f"""
    Error when merging partial logs for the region. Please
    check that pre-existing results are consistent with each other."""

    # Keys that must be identical
    identical_keys = [
        "log10p_threshold", "regenie_columns", "csv_columns",
        "synapse_folder_id"
    ]

    # Keys for lists that must have same content (maybe different order)
    equivalent_keys = ["all_tar_files"]

    # Those keys become non-empty only after merging
    empty_list_keys = [
        "all_proteins", "kept_proteins", "partial_log_merged",
        "partial_output_filenames"
    ]

    # Keys that must be the union but that must contain no duplicates
    union_no_dups_keys = [
        "tar_processed", "tar_downloaded_skip",
    ]

    # Keys that may have duplicates
    union_keys = ["tar_skipped"]

    # Keys that must be different between each log and not None
    different_keys = ["log_filename", "output_filename"]

    # Keys to treat differently but that are not protein names
    other_keys = []

    # Non-protein keys
    non_protein_keys = [
        *identical_keys, *equivalent_keys, *empty_list_keys,
        *union_no_dups_keys, *union_keys, *different_keys, *other_keys
    ]

    # ----------- Keys that must be identical -----------
    if len(compatible_log_files) > 0:
        first_log_dict = next(iter(compatible_log_files.values()))

        for key in identical_keys:

            assert all([
                log_dict[key] == first_log_dict[key]
                for log_dict in compatible_log_files.values()
            ])

            # Take the first log as reference
            merged_log[key] = first_log_dict[key]

    # ----------- Keys that must be equivalent -----------
    if len(compatible_log_files) > 0:
        first_log_dict = next(iter(compatible_log_files.values()))

        for key in equivalent_keys:

            # Same content (same unique elements)
            assert all([
                set(log_dict[key]) == set(first_log_dict[key])
                for log_dict in compatible_log_files.values()
            ])
            # Same number of elements
            assert all([
                len(log_dict[key]) == len(first_log_dict[key])
                for log_dict in compatible_log_files.values()
            ])

            # Take the first log as reference
            merged_log[key] = first_log_dict[key]
            merged_log[key].sort()

    # ----------- Keys that must be different -----------
    for key in different_keys:

        # Check that there is no duplicate if we take the list of all values
        all_values = [
            log_dict[key]
            for log_fname, log_dict in compatible_log_files.items()
        ]
        msg_err = f"[ERR] {key} should be different in all log files, but found duplicates.{merging_msg}"
        assert len(set(all_values)) == len(all_values), f"{msg_err}"

        # We don't accept None, all must be strings
        msg_err = f"[ERR] {key} can not contain None values.{merging_msg}"
        assert None not in all_values, f"{msg_err}"

        # Then we leave values untouched

    # ----- Keys that should be the union (but with no duplicates) -----
    for key in union_no_dups_keys:

        all_values = []
        for log_fname, log_dict in compatible_log_files.items():
            all_values += log_dict[key]

        # Check that there is no duplicate if we take the list of all values
        n_tot = len(all_values)
        all_values = list(set(all_values))
        n_unique = len(all_values)
        msg_err =  f"[ERR] Duplicates found for key {key}"
        assert n_unique == n_tot, f"{msg_err}"

        # Then we take the union of all values
        merged_log[key] = all_values
        merged_log[key].sort()

    # ----- Keys that should be the union (may have duplicates) -----
    for key in union_keys:

        all_values = []
        for log_fname, log_dict in compatible_log_files.items():
            all_values += log_dict[key]

        # There could be duplicates, but we don't keep only unique values
        all_values = list(set(all_values))

        merged_log[key] = all_values
        merged_log[key].sort()

    # --------------------- tar skipped --------------------------------
    # The true 'tar_skipped' remaining is a subset of the union of all
    # tar_skipped, minus the union of tar_processed
    merged_log["tar_skipped"] = [
        tar_fname for tar_fname in merged_log["tar_skipped"]
        if tar_fname not in merged_log["tar_processed"]
    ]

    # --------------------- Protein keys -------------------------------
    # Then merged all the remaining keys, which should be protein names
    for log_fname, log_dict in compatible_log_files.items():

        # Protein in this log file
        remaining_keys = [
            k for k in log_dict.keys()
            if k not in non_protein_keys
        ]
        # Sort them so that we write them in json sorted
        remaining_keys.sort()

        for key in remaining_keys:
            msg_err =  f"[ERR] Protein {key} found in multiple log files.{merging_msg}"
            assert key not in merged_log, f"{msg_err}"
            merged_log[key] = log_dict[key]

    # ------------------ Merged log files references -------------------

    # Update the list of output and log fnames
    concat_keys = [
        ("partial_log_merged", "log_filename"),
        ("partial_output_filenames","output_filename")
    ]

    for key_merged, key_source in concat_keys:
        merged_log[key_merged] = []
        for log_fname, log_dict in compatible_log_files.items():
            merged_log[key_merged].append(log_dict[key_source])

        merged_log[key_merged].sort()

        # Check that there is no dup
        n_tot = len(merged_log[key_merged])
        merged_log[key_merged] = list(set(merged_log[key_merged]))
        n_unique = len(merged_log[key_merged])
        msg_err = f"[ERR] Duplicates found in {key_source} between logs. {merging_msg}"
        assert n_unique == n_tot, f"{msg_err}"


    return merged_log

def merge_partial_output_files(
    output_fname: str,
    partial_output_filenames: list[str],
) -> None:
    """
    Concatenate partial text output files into a single output file.

    This function will create a new text file.

    Parameters
    ----------
    output_fname : str
        Path to the final merged output text file to create, e.g.
        ``"syn51365308-output_text-2026-05-12--10:00:00.txt"``.
    partial_output_filenames : list[str]
        Paths to the partial output files to concatenate, e.g.
        ``["PART-syn51365308-output_text-2026-05-12--09:00:00.txt",
        ...]``. Files are sorted by filename before concatenation so
        that order matches creation order. Each file's content is
        wrapped with a separator banner showing the source filename.
    """

    # We want them sorted because it matches the creation order
    partial_output_filenames.sort()

    with open(output_fname, 'wt') as outfile:
        for fname in partial_output_filenames:
            with open(fname) as infile:

                # Add a separator between each part file
                outfile.write("─"*80 + "\n")
                outfile.write(f"{' '*6} From: {os.path.basename(fname)} " + "\n")
                outfile.write("─"*80 + "\n\n")

                # Add content of part file
                outfile.write(infile.read())

        # And a final separator
        outfile.write("\n\n" + "─"*80 + "\n")
        outfile.write(f"{' '*6} From: {os.path.basename(output_fname)} " + "\n")
        outfile.write("─"*80 + "\n\n")



def keep_significant_qtls_from_region(
        synapse_folder_id: str = PGWAS_REGIONS["Combined"],
        download_location: str = "./ukb_ppp_dl/data",
        res_location: str = "./ukb_ppp_dl/results",
        login_kwargs: dict[str, Any] = {},
        regenie_sep: str = " ",
        regenie_columns: list[str] | None = None,
        csv_columns: list[str] | None = None,
        log10p_threshold: float | int = 7,
        create_log: bool | int = True,
        log_kwargs: dict[str, Any] = {},  # Not fully implemented yet,
        protein_to_process: list[str] | None = None,
        verbose: bool | int = False,
        delete_downloaded_tar: bool = True,
        delete_chr_csv: bool = True,
        delete_tar_csv: bool = True,
        delete_tar_log: bool = True,
        delete_partial_logs: str | bool = "current",
        delete_partial_outputs: str | bool = "current",
) -> tuple[pl.DataFrame, dict]:
    """
    Keep significant pGWAS QTLs for an entire Synapse region folder.

    This is the main entry point for processing one ancestry group. It:

    1. Lists all protein tar files in the Synapse folder.
    2. For each protein, downloads the tar file from Synapse (if not
       already downloaded), extracts per-chromosome REGENIE files,
       filters for significant QTLs, and merges them into a
       protein-level CSV.
    3. Concatenates all protein-level results into one region-level CSV.
    4. Documents the process and results in logs and output text files.
    5. Optionally cleans up intermediate files.

    Partial results from a previous interrupted run are automatically
    detected and reused.

    This function does not require a lot of free space to run because it
    processes one protein at a time and deletes intermediate files along
    the way (if specified).

    This function creates logs and output text files for each run. This
    is important and has several benefits:

    1. It allows the function to be safely re-run multiple times without
       overwriting previous logs and outputs.
    2. It allows the function to automatically detect and reuse partial
       runs
    3. It provides a detailed record of the parameters used to create
       the generated results. This can be crucial for reproducibility in
       science.
    4. It facilitates sanity checks and monitoring of the results.

    It is highly recommended to keep the logs and output text files (at
    least at the region level, that is to say with ``create_log=True``
    or ``create_log=1``) for future reference.

    Note that even if ``create_log=False``, partial log files will still
    be created, but they will be deleted at the end of the function.

    This function will create many files during the process, but most of
    them are intermediate files that can be deleted at the end of the
    run if specified. The main final output is the region-level CSV file
    containing all significant QTLs across all processed proteins, and
    the region-level output file documenting the process and results.

    Parameters
    ----------
    synapse_folder_id : str, optional
        Synapse folder ID for the GWAS region, e.g. ``"syn51365308"``
        for the Combined ancestry. Defaults to
        ``PGWAS_REGIONS["Combined"]``.
    download_location : str, optional
        Local directory where tar files are downloaded.
    res_location : str, optional
        Local directory where result CSV, log and output files are
        written.
    login_kwargs : Dict[str, Any], optional
        Keyword arguments forwarded to
        ``synapseclient.Synapse.login()``. Common keys include
        ``"authToken"`` (``str``), ``"email"`` (``str``) and
        ``"profile"`` (``str``). See `synapseclient.Synapse.login()
        <https://python-docs.synapse.org/en/stable/reference/client/#synapseclient.Synapse.login>`_
        for more information.

        Alternatively, you can configure your Synapse credentials in a
        ``.synapseConfig`` file and leave this argument empty to use the
        default login behaviour. You can find more information about
        ``.synapseConfig`` in the `Synapse documentation
        <https://python-docs.synapse.org/en/stable/tutorials/authentication/>`_.
    regenie_sep : str, optional
        Field separator used in REGENIE files. Default is ``" "``
        (space).
    regenie_columns : list[str] or None, optional
        Subset of REGENIE columns to load, e.g. ``["CHROM", "GENPOS",
        "ID", "BETA", "SE", "LOG10P"]``. ``"LOG10P"`` must be included.
        ``None`` loads all columns: ``["CHROM", "GENPOS", "ID",
        "ALLELE0", "ALLELE1", "A1FREQ", "INFO", "N", "TEST", "BETA",
        "SE", "CHISQ", "LOG10P", "EXTRA"]``.
    csv_columns : list[str] or None, optional
        New column names to assign after loading. Must match the length
        of *regenie_columns* when both are provided. ``"LOG10P"`` must
        keep the same name. If ``None``, no renaming is done and
        original REGENIE column names are used.
    log10p_threshold : float or int, optional
        Minimum ``LOG10P`` value for significance. Default is ``7``.
    create_log : bool or int, optional
        Logging verbosity. ``0`` disables; ``1`` creates region-level
        logs; ``2`` also creates protein-level logs; ``3`` also creates
        per-chromosome logs. Default is ``True`` (``1``).
    log_kwargs : dict[str, Any], optional
        Keyword arguments forwarded to :func:`save_log`. Recognised key:
        ``"log_fname"`` (``str``) to override the default log filename.
        This parameter is not fully implemented yet and it is mostly the
        default filenames that are used.
    protein_to_process : list[str] or None, optional
        Whitelist of proteins to process. Accepts Synapse IDs (e.g.
        ``"syn52361344"``) or tar filenames (e.g.
        ``"ABCA2_Q9BZC7_OID30146_v1_Cardiometabolic_II.tar"``). ``None``
        processes all proteins.
    verbose : bool or int, optional
        Verbosity level. ``0``/``False`` is silent; higher values print
        more details at each processing step.
    delete_downloaded_tar : bool, optional
        Delete each downloaded tar file after processing to save disk
        space. Default is ``True``.
    delete_chr_csv : bool, optional
        Delete per-chromosome significant-QTL CSV files after merging
        them into the protein-level CSV. Default is ``True``.
    delete_tar_csv : bool, optional
        Delete protein-level significant-QTL CSV files after merging
        them into the region-level CSV. Default is ``True``.
    delete_tar_log : bool, optional
        Delete protein-level log files after the region-level log is
        created. Default is ``True``.
    delete_partial_logs : str or bool, optional
        Controls deletion of partial region log files. ``"current"`` or
        ``True`` deletes only the log created by this run; ``"all"``
        deletes all compatible partial logs; ``False`` keeps them all.
        Default is ``"current"``.
    delete_partial_outputs : str or bool, optional
        Controls deletion of partial region output text files. Same
        options as *delete_partial_logs*. Default is ``"current"``.

    Returns
    -------
    all_significant_qtls : pl.DataFrame
        DataFrame of all significant QTLs across all processed proteins.
        Contains a leading ``"protein_name"`` (``str``) column followed
        by the REGENIE (or renamed) columns.
    final_log_reg : dict
        Final region log dictionary (merged from all partial logs).
    """

    # NOTE: In any case, we are producing an output text file
    full_date = datetime.today().strftime('%Y-%m-%d--%H:%M:%S')
    fname_base = f'{res_location}/{synapse_folder_id}-significant_qtls-{full_date}'
    fname_base_PART = f'{res_location}/PART-{synapse_folder_id}-significant_qtls-{full_date}'

    # Make sure path exists otherwise create it
    p = pathlib.Path(f"{fname_base}-text.txt")
    p.parent.mkdir(parents=True, exist_ok=True)

    fout = open(f"{fname_base_PART}-text.txt", 'wt')
    sys.stdout = fout

    t_start_region = time.time()

    # ----------- List of protein tar files in Region folder -----------
    # Get a list of (synapse_id, tar_name) for all the protein tar files
    folder_tar_entities = list_tar_files_in_region_folder(
        synapse_id=synapse_folder_id,
        login_kwargs=login_kwargs,
    )

    # Keeping only a subset if protein_to_process is specified
    tar_entities = [
        (synapse_id, tar_name)
        for synapse_id, tar_name in folder_tar_entities
        if (
            (protein_to_process is None) or
            (synapse_id in protein_to_process) or
            (tar_name in protein_to_process)
        )
    ]

    if int(verbose) > 0:
        print(f"Found {len(folder_tar_entities)} protein tar files in Synapse folder {synapse_folder_id}.")
        print(f"Keeping only {len(tar_entities)} protein tar files after filtering for specified proteins.", flush=True)

    # ----------------- Preparing (partial) log for the region ---------

    dict_df_significant = {}
    dict_csv_significant = {}

    log_reg = {
        "synapse_folder_id": synapse_folder_id,
        "log10p_threshold": log10p_threshold,
        "regenie_columns": regenie_columns,
        "csv_columns": csv_columns,
        "all_tar_files": [tar_name for _, tar_name in tar_entities],
        "log_filename": None,
        "output_filename": f"{fname_base_PART}-text.txt",
        "tar_skipped" : [],
        "tar_processed": [],
        "tar_downloaded_skip": [],
        "all_proteins": [],
        "kept_proteins" : [],
        "partial_log_merged": [],
        "partial_output_filenames" : [],
    }

    # ---- Find pre-existing partial region log files ------------------
    compatible_log_files = find_partial_region_logs(
        example_log_dict = log_reg,
        res_location=res_location,
        verbose=int(verbose) - 1,
    )

    # This dict can be empty (so expected keys are not present, be careful!))
    merged_partial_region_log = merge_partial_region_logs(compatible_log_files)

    # ---- Save current partial region log files ------------------

    # Save initial part log
    save_log(
        f"{fname_base_PART}-log.json", log_reg,
        overwrite=False, add_date=False, new_name=True, verbose=False,
    )

    # --------------  Process each protein one by one -----------------
    for synapse_id, tar_name in tar_entities:

        # Preparing result filename keeping only significant QTLs
        protein_name = f"{Path(tar_name).stem.split('_')[0]}"
        res_merged_fname = f"{res_location}/{protein_name}-significant_qtls"

        if int(verbose) > 0:
            print(f"\n┌─{'─'*62}─┐")
            print(f"{" "*10} Protein {protein_name} | Synapse ID: {synapse_id}")
            print(f"\n└─{'─'*62}─┘", flush=True)


        # -------- Skipping if results pre-exist --------------

        # Don't run again if result file already exists
        # and if a corresponding partial log file corresponds
        # Just load the results.
        if (
            os.path.isfile(f"{res_merged_fname}.csv") and
            tar_name in merged_partial_region_log.get("tar_processed", [])
        ):
            log_reg["tar_skipped"].append((tar_name))
            dict_df_significant[protein_name] = pl.read_csv(f"{res_merged_fname}.csv")
            dict_csv_significant[protein_name] = f"{res_merged_fname}.csv"

            if int(verbose) > 0:
                print(
                    f"Result file already created: {res_merged_fname}.csv. Loading results.",
                    flush=True
                )

            if int(create_log) > 0:
                if int(verbose) > 0:
                    print(f"Updating partial log file for the region: {f"{fname_base_PART}-log.json"}")
                save_log(
                    f"{fname_base_PART}-log.json", log_reg,
                    overwrite=True, add_date=False
                )

        else:

            t_start = time.time()

            # ------------ Download tar file --------------

            tar_local_fname, skipped = download_from_synapse(
                synapse_id= synapse_id,
                download_location=download_location,
                login_kwargs=login_kwargs,
                verbose=verbose,
            )

            if skipped:
                log_reg["tar_downloaded_skip"].append(tar_name)

            # ------------ Process tar file --------------

            all_csv_fnames, log_tar = process_one_tar_file(
                tar_local_fname,
                res_location=res_location,
                separator=regenie_sep,
                columns = regenie_columns,
                new_columns=csv_columns,
                log10p_threshold=log10p_threshold,
                create_log=int(create_log) - 1,
                log_kwargs=log_kwargs,
                verbose=int(verbose) - 1,
            )

            # --- Merge chr csv files that kept only significant QTLs --

            df_significant_qtls_prot, log_merged = merge_significant_qtls_from_csv(
                all_csv_fnames,
                output_fname=f"{res_merged_fname}.csv",
                create_log=False, # create log after merge
                log_kwargs=log_kwargs,
                delete_csv=delete_chr_csv,
                verbose=int(verbose) - 1,
            )

            dict_df_significant[protein_name] = df_significant_qtls_prot
            dict_csv_significant[protein_name] = f"{res_merged_fname}.csv"

            log_reg["tar_processed"].append(tar_name)

            t_end = time.time()

            if int(verbose) > 0:
                print(f"Processed tar file {tar_name} in {t_end - t_start:.2f} s", flush=True)

            # ----------- Update logs for protein -----------
            # Merge dict (Duplicate keys like log_filename: keep log_tar value)
            log_tar = log_merged | log_tar
            log_tar["synapse_id"] = synapse_id
            log_tar["time_taken"] = t_end - t_start


            if int(create_log) > 1:
                # --------- Create log file for the protein ---------------
                log_fname = f"{res_merged_fname}-log.json"

                # There might be already a log file from the call to
                # process_one_tar_file function, we will delete it
                # after creating the new tar log file with more info
                old_log_tar_fname = log_tar.pop("log_filename", None)

                save_log(log_fname, log_tar, **log_kwargs)

                # Delete old tar log file
                if old_log_tar_fname is not None and os.path.isfile(old_log_tar_fname):
                    os.remove(old_log_tar_fname)

            #  -------- Update partial log file for the region ----------

            # Don't include keys that are common to all tar file
            # And don't include keys that are "too much"
            log_reg[protein_name] = {
                k: v for k, v in log_tar.items()
                if k not in ["skipped_chr_files", "all_csv_fnames", "log10p_threshold", "regenie_columns", "csv_columns"]
            }
            log_reg[protein_name]["n_skipped_chr_files"] = len(log_tar["skipped_chr_files"])
            log_reg[protein_name]["time_taken"] = t_end - t_start

            if int(create_log) > 0:
                if int(verbose) > 0:
                    print(f"Updating partial log file for the region: {f"{fname_base_PART}-log.json"}")
                save_log(
                    f"{fname_base_PART}-log.json", log_reg,
                    overwrite=True, add_date=False, verbose=False
                )


            # ----- Delete the downloaded tar file to save space (optional)
            if not skipped and delete_downloaded_tar:
                if int(verbose) > 0:
                    print(f"Deleting downloaded tar file {tar_local_fname} to save space.")
                os.remove(tar_local_fname)

    # ----------------- Verbose ---------------------
    if log_reg["tar_skipped"] and int(verbose) > 0:
        print(f"[WARNING]: {len(log_reg['tar_skipped'])}/{len(log_reg['all_tar_files'])}  skipped tar files! ")
        print(_CONSISTENCY_MESSAGE)

    # ------- Add one column to all df with the protein name -----------

    # Add one column with the protein name to each protein df
    for protein in dict_df_significant.keys():

        # To be able to reorder columns to have protein_name afterwards
        original_columns = dict_df_significant[protein].columns

        # Add the column
        dict_df_significant[protein] = dict_df_significant[protein].with_columns(
            protein_name = pl.lit(protein)
        )

        # Reorder columns to have protein_name first
        dict_df_significant[protein] = dict_df_significant[protein].select(
            ["protein_name"] + original_columns
        )

    # ------------ Merge partial logs with current log -------

    # Merge the current partial log with the merged pre-existing region logs
    compatible_log_files[log_reg["log_filename"]] = log_reg

    final_log_reg = merge_partial_region_logs(compatible_log_files)

    # ------------- Merge all merged protein dataframes ----------------

    if verbose:
        print(f"{"="*80}")
        print(f"Merging significant QTLs from {len(dict_df_significant)} protein files.", flush=True)

    # Distinguish between proteins with and without at least one significant QTL
    final_log_reg["kept_proteins"] = [
        prot for prot, df in dict_df_significant.items() if len(df) > 0
    ]
    final_log_reg["all_proteins"] = list(dict_df_significant.keys())

    # Concat all the protein dataframes
    non_empty_dfs = [
        df for prot, df in dict_df_significant.items() if len(df) > 0
    ]

    all_significant_qtls = pl.concat(non_empty_dfs)

    # Save the concatenated dataframe to a new file
    all_significant_qtls.write_csv( f"{fname_base}.csv")

    if verbose:
        print(f"Found {len(final_log_reg['kept_proteins'])}/{len(final_log_reg['all_proteins'])} proteins with at least one significant QTL.")
        print(f"All significant QTLs saved to:  {f"{fname_base}.csv"}")
        print(f"Total significant QTLs:         {len(all_significant_qtls)}")
        print(f"{"="*80}", flush=True)

    # ------------------- Merge output files ---------------------------
    fout.close()

    # Merge all partial output files into "output_fname" file
    merge_partial_output_files(
        f"{fname_base}-text.txt", final_log_reg["partial_output_filenames"]
    )

    # Use this final output file as the output filename (in append mode!)
    fout = open(f"{fname_base}-text.txt", 'at')
    sys.stdout = fout
    final_log_reg["output_filename"] = f"{fname_base}-text.txt"

    # ---------------Create log file for the region ------------------
    if int(create_log) > 0:
        # If log_fname not in log_kwargs, use a default one
        save_log(
            f"{fname_base}-log.json", final_log_reg,
            overwrite=False, add_date=False, new_name=True, verbose=False,
        )

    # ---------------------Sanity Checks ------------------

    # At the end, there should be no skipped tar file
    err_msg = f"[ERR] There are still skipped tar files after merging logs: {final_log_reg["tar_skipped"]})."
    assert final_log_reg["tar_skipped"] == [], err_msg

    # At the end, all tar files should be in the list of processed tar files
    err_msg = f"[ERR] Missing tar files in the list of processed tar files."
    expected_tar_files = set(tar_name for tar_id, tar_name in tar_entities)
    assert set(final_log_reg["tar_processed"]) == expected_tar_files, f"{err_msg}. Expected: {expected_tar_files}, got: {set(final_log_reg['tar_processed'])}"
    assert set(final_log_reg["tar_processed"]) == set(final_log_reg["all_tar_files"]), f"{err_msg}. Expected: {set(final_log_reg['all_tar_files'])}, got: {set(final_log_reg['tar_processed'])}"

    # At the end, the number of proteins and tar files must be consistent
    err_msg = f"[ERR] The number of processed tar files does not correspond to the number of proteins in the final log."
    assert len(final_log_reg["tar_processed"]) == len(final_log_reg["all_proteins"]), f"{err_msg} Expected: {len(final_log_reg['tar_processed'])}, got: {len(final_log_reg['all_proteins'])}"

    # ---------------- Delete tar csv files ----------------------
    if delete_tar_csv:
        fnames = [
            final_log_reg[protein]["merged_csv_filename"]
            for protein in final_log_reg["all_proteins"]
        ]

        n_deleted = delete_files(fnames, verbose=int(verbose) - 1)
        if n_deleted > 0 and int(verbose) > 0:
            print(f"Deleted {n_deleted} protein CSV files.")

    # ---------------- Delete tar log files ----------------------
    if delete_tar_log:
        fnames = [
            final_log_reg[protein]["log_filename"]
            for protein in final_log_reg["all_proteins"]
        ]
        n_deleted = delete_files(fnames, verbose=int(verbose) - 1)
        if n_deleted > 0 and int(verbose) > 0:
            print(f"Deleted {n_deleted} protein log files.")

    # ---------------- Delete partial log files ----------------------
    # Deleting all partial logs used to create this final log
    if delete_partial_logs == "all":
        fnames = final_log_reg["partial_log_merged"]
        n_deleted = delete_files(fnames, verbose=int(verbose) - 1)
        if int(verbose) > 0 and n_deleted > 0:
            print(f"Deleted {n_deleted} partial log files.")

    # Deleting only the current partial log created by this run
    elif delete_partial_logs == "current" or delete_partial_logs == True:
        if os.path.isfile(f"{fname_base_PART}-log.json"):
            os.remove(f"{fname_base_PART}-log.json")
            if int(verbose) > 0:
                print(f"Deleted current partial log file {f'{fname_base_PART}-log.json'}.")

    # ---------------- Delete partial output files ----------------------
    # Deleting all partial outputs used to create this final output
    if delete_partial_outputs == "all":
        fnames = final_log_reg["partial_output_filenames"]
        n_deleted = delete_files(fnames, verbose=int(verbose) - 1)
        if int(verbose) > 0 and n_deleted > 0:
            print(f"Deleted {n_deleted} partial output files.")

    # Deleting only the current partial output created by this run
    elif delete_partial_outputs == "current" or delete_partial_outputs == True:
        if os.path.isfile(f"{fname_base_PART}-text.txt"):
            os.remove(f"{fname_base_PART}-text.txt")
            if int(verbose) > 0:
                print(f"Deleted current partial output file {f'{fname_base_PART}-text.txt'}.")

    # ------------------- Verbose -----------------------------------
    t_end_region = time.time()
    if int(verbose) > 0:
        print(f"\n{"="*80}", flush=True)
        print(f"Processed region tar files in {t_end_region - t_start_region:.2f} s", flush=True)

    if final_log_reg["tar_downloaded_skip"] and int(verbose) > 0:
        print(f"{len(final_log_reg['tar_downloaded_skip'])}/{len(final_log_reg['all_tar_files'])} tar files were already downloaded.")


    fout.close()
    return all_significant_qtls, final_log_reg