"""
General design choices

- The code checks if the result file already exists before running the analysis, to avoid re-running it if it has already been run before.
- If you don't want to reuse result files, you must delete them yourself manually.
"""
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

from .common import save_log, delete_files, download_from_synapse

CONSISTENCY_MESSAGE = """
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
        synapse_id=PGWAS_REGIONS["Combined"],
        login_kwargs={},
    ):
    """
    List available protein tar files in a Synapse directory.

    Parameters:
    synapse_id (str): The Synapse ID of the directory to list files from.

    Returns:
    list: A list of Synapse IDs and names of the available protein tar files.
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
        chr_file_gz,
        separator=" ",
        columns = None,
        new_columns=None,
        log10p_threshold=7,
        create_log = False,
        log_kwargs = {},
        verbose=False,
    ):
    # NOTE: LOG10P column must exist in the dataframe and with the correct case

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
        protein_tar_file,
        chr_gz_fname,
        res_location="./ukb_ppp_dl/results",
        separator=" ",
        columns = None,
        new_columns=None,
        log10p_threshold=7,
        create_log = False,
        log_kwargs = {},
        verbose=False,
    ):

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
        csv_fnames,
        output_fname = None,
        create_log = False,
        log_kwargs = {},
        delete_csv = False,
        verbose = False,
    ):
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
        tar_fname,
        res_location="./ukb_ppp_dl/results",
        separator=" ",
        columns = None,
        new_columns=None,
        log10p_threshold=7,
        create_log = True,
        log_kwargs = {},
        verbose = False,
    ):
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
        print(CONSISTENCY_MESSAGE)

    if int(create_log) > 0:
        # If log_fname not in log_kwargs, use a default one
        log_fname = f"{res_location}/{Path(tar_fname).stem}-log.json"
        log_fname = log_kwargs.pop("log_fname", log_fname)
        save_log(log_fname, log_tar, **log_kwargs)

    return log_tar["all_csv_fnames"], log_tar

def find_partial_region_logs(
    example_log_dict = None,
    synapse_folder_id=PGWAS_REGIONS["Combined"],
    res_location='./ukb_ppp_dl/results',
    regenie_columns = None,
    csv_columns = None,
    log10p_threshold = 7,
    all_tar_files = None,
    verbose = False,
):

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

        if all_tar_files is not None and set(log_dict["all_tar_files"]) != set(all_tar_files):
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
    compatible_log_files,
):
    # This function should be called only for log files created with the same parameters and created by running the function process_one_region_folder
    # We cannot use this function with already merged log files.


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
    output_fname,
    partial_output_filenames,
) :

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



def process_one_region_folder(
        synapse_folder_id=PGWAS_REGIONS["Combined"],
        download_location="./ukb_ppp_dl/data",
        res_location="./ukb_ppp_dl/results",
        login_kwargs={},
        regenie_sep = " ",
        regenie_columns = None,
        csv_columns = None,
        log10p_threshold = 7,
        create_log = True,
        log_kwargs = {},
        protein_to_process = None,
        verbose=False,
        delete_downloaded_tar = True,
        delete_chr_csv = True,
        delete_tar_csv = True,
        delete_tar_log = True,
        delete_partial_logs = "current",
        delete_partial_outputs = "current",
):
    # NOTE: Even if create_log is False, partial log files will be created
    # but they will be deleted at the end of the function if create_log is False

    # NOTE: In any case, we are producing an output text file

    full_date = datetime.today().strftime('%Y-%m-%d--%H:%M:%S')
    output_fname = f'{res_location}/{synapse_folder_id}-output_text-{full_date}.txt'
    part_output_fname = f'{res_location}/PART-{synapse_folder_id}-output_text-{full_date}.txt'
    res_fname = f"{res_location}/{synapse_folder_id}-significant_qtls"
    part_resfname = f"{res_location}/PART-{synapse_folder_id}-significant_qtls"

    # Make sure path exists otherwise create it
    p = pathlib.Path(output_fname)
    p.parent.mkdir(parents=True, exist_ok=True)

    fout = open(part_output_fname, 'wt')
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
        "output_filename": part_output_fname,
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

    part_log_reg_fname = f"{part_resfname}-log-{full_date}.json"
    part_log_reg_fname = log_kwargs.pop("log_fname", part_log_reg_fname)

    # Save initial part log
    save_log(
        part_log_reg_fname, log_reg,
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
                    print(f"Updating partial log file for the region: {part_log_reg_fname}")
                save_log(
                    part_log_reg_fname, log_reg,
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
                    print(f"Updating partial log file for the region: {part_log_reg_fname}")
                save_log(
                    part_log_reg_fname, log_reg,
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
        print(CONSISTENCY_MESSAGE)

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
    all_significant_qtls.write_csv( f"{res_fname}.csv")

    if verbose:
        print(f"Found {len(final_log_reg['kept_proteins'])}/{len(final_log_reg['all_proteins'])} proteins with at least one significant QTL.")
        print(f"All significant QTLs saved to:  {f"{res_fname}.csv"}")
        print(f"Total significant QTLs:         {len(all_significant_qtls)}")
        print(f"{"="*80}", flush=True)


    # ------------------- Merge output files ---------------------------
    fout.close()

    # Merge all partial output files into "output_fname" file
    merge_partial_output_files(
        output_fname, final_log_reg["partial_output_filenames"]
    )

    # Use this final output file as the output filename (in append mode!)
    fout = open(output_fname, 'at')
    sys.stdout = fout
    final_log_reg["output_filename"] = output_fname

    # ---------------Create log file for the region ------------------
    if int(create_log) > 0:
        # If log_fname not in log_kwargs, use a default one
        save_log(f"{res_fname}-log.json", final_log_reg, **log_kwargs)


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
        if os.path.isfile(part_log_reg_fname):
            os.remove(part_log_reg_fname)
            if int(verbose) > 0:
                print(f"Deleted current partial log file {part_log_reg_fname}.")

    # ---------------- Delete partial output files ----------------------
    # Deleting all partial outputs used to create this final output
    if delete_partial_outputs == "all":
        fnames = final_log_reg["partial_output_filenames"]
        n_deleted = delete_files(fnames, verbose=int(verbose) - 1)
        if int(verbose) > 0 and n_deleted > 0:
            print(f"Deleted {n_deleted} partial output files.")

    # Deleting only the current partial output created by this run
    elif delete_partial_outputs == "current" or delete_partial_outputs == True:
        if os.path.isfile(part_output_fname):
            os.remove(part_output_fname)
            if int(verbose) > 0:
                print(f"Deleted current partial output file {part_output_fname}.")

    # ------------------- Verbose -----------------------------------
    t_end_region = time.time()
    if int(verbose) > 0:
        print(f"\n{"="*80}", flush=True)
        print(f"Processed region tar files in {t_end_region - t_start_region:.2f} s", flush=True)

    if final_log_reg["tar_downloaded_skip"] and int(verbose) > 0:
        print(f"{len(final_log_reg['tar_downloaded_skip'])}/{len(final_log_reg['all_tar_files'])} tar files were already downloaded.")


    fout.close()
    return all_significant_qtls, log_reg