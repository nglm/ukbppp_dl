"""
General design choices

- The code checks if the result file already exists before running the analysis, to avoid re-running it if it has already been run before.
- If you don't want to reuse result files, you must delete them yourself manually.
- The log files stores information about the results that won't be accessible outside the scope of a given function.
- The log file stores information about parameters used only when the parameteres are not directly used when calling the function. In other words, if a user calls a given function, the parameters used in the call won't be saved i n the log file.
"""

import time
import synapseclient
from synapseclient.models import Folder
import tarfile
import gzip
import polars as pl
import os
import json
import pathlib

# Synapse directory containing pQTL summary statistics (here for Europe)
REGION_PQTL_DIR = 'syn51365303'

DOWNLOAD_LOCATION = "./data/"

# Specific protein tar files for test purpose
ACOT13_FILE = "syn52362654"
ZNF174_FILE = "syn52363271"

# Mandatory columns in the .regenie files
MANDATORY_COLUMNS = ["CHROM", "GENPOS", "ID", "BETA", "SE", "LOG10P"]
NEW_COLUMN_NAMES = ["chrom", "qtl_pos", "qtl_id", "beta", "se", "log10p"]

# Significance threshold for pQTLs
LOG10P_THRESHOLD = 7

# Separator used in the .regenie files
REGENIE_SEP = " "

# Whether to create a log file
CREATE_LOG = False

# Synapse login kwargs
LOGIN_KWARGS = {}

CONSISTENCY_MESSAGE = """
    Files were skipped because their result file already exists.
    Please make sure that the pre-existing result files are correct and consistent with the new analysis.
    If not consistent, delete result files and re-run the analysis.
    Decisive parameters include: log10p_threshold, regenie_columns, csv_columns.
"""


def create_log(log_fname, log_dict, overwrite=True, verbose=False):

    if os.path.isfile(log_fname) and not overwrite and verbose:
        print(f"Log file {log_fname} already exists. Not overwriting.")
        return False

    else:

        # Make sure path exists otherwise create it
        p = pathlib.Path(log_fname)
        p.parent.mkdir(parents=True, exist_ok=True)

        with open(log_fname, 'w') as f_log:
            json.dump(log_dict, f_log, indent=2)

        if verbose:
            print(f"Saved log file to {log_fname}.")

        return True

def list_available_protein_tar_files(
        synapse_id=REGION_PQTL_DIR,
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

    # This will only get immediate children, not subfolders' contents
    folder = Folder(id=synapse_id)
    folder = folder.sync_from_synapse(download_file=False, recursive=False)

    # Filter for tar files and return their Synapse IDs and names
    tar_files = [
        (file.id, file.name)
        for file in folder.files if file.name.endswith('.tar')
    ]

    return tar_files

def download_protein_tar_file(
        synapse_id,
        download_location=".",
        login_kwargs={},
        verbose=False,
    ):
    """
    Download a protein tar file from Synapse.

    If download_location is not specified, the file will be downloaded to the synapse cache (see ~/.synapseConfig or ~/.synapseCache)

    Parameters:
    synapse_id (str): The Synapse ID of the file to download.
    download_location (str): The directory where the file will be downloaded. Default is the current directory.

    Returns:
    str: The path to the downloaded file.
    """
    syn = synapseclient.Synapse()
    # Use your synapse account here, as configured in your .synapseConfig file.
    syn.login(**login_kwargs)

    # ----- First check if the file has already been downloaded ------
    # If downloadFile is False, the file will not be downloaded, but we have access to synapse.entity.Entity attributes such as name.
    entity = syn.get(synapse_id, downloadFile=False)

    # Check if the result file already exists in the download location
    expected_file_path = f"{download_location}/{entity.name}"
    if os.path.isfile(expected_file_path):
        skipped = True
        if verbose:
            print(f"File already downloaded: {expected_file_path}. Skipping download.")
    else:
        skipped = False
        # Download the file to the specified location
        entity = syn.get(
            synapse_id, downloadLocation=download_location, downloadFile=True
        )
        expected_file_path = entity.path
        if verbose:
            print(f"Downloaded file: {expected_file_path}")

    return expected_file_path, skipped

def keep_significant_qtls_from_chr_gz_file(
        chr_file_gz,
        separator=REGENIE_SEP,
        columns = MANDATORY_COLUMNS,
        new_columns=NEW_COLUMN_NAMES,
        log10p_threshold=LOG10P_THRESHOLD,
        verbose=False,
    ):

    # Opening the extracted .gz file to read its content
    with gzip.open(chr_file_gz, 'rt') as f_regenie:

        # Reading the .regenie file as a CSV and creating a dataframe
        summary_stats = pl.read_csv(
            f_regenie,
            separator=separator,
            columns=columns,
            new_columns=new_columns,
            )


        n_tot_qtls = pl.len(summary_stats)
        all_qtls = summary_stats.get_column("qtl_id").to_list()

        # Filtering the summary statistics to keep only significant QTLs
        summary_stats_significant = summary_stats.filter(
            pl.col("log10p") >= log10p_threshold
        )

        n_kept_qtls = pl.len(summary_stats_significant)

        log = {
            "n_tot_qtls": n_tot_qtls,
            "n_kept_qtls": n_kept_qtls,
            "all_qtls": all_qtls,
        }

        if verbose:
            print(f"Total QTLs in chr file: {n_tot_qtls}")
            print(f"Kept QTLs (>= {log10p_threshold}):       {n_kept_qtls}")

    return summary_stats_significant, log

def process_one_chr_from_protein_tar_file(
        protein_tar_file,
        chr_gz_fname,
        separator=REGENIE_SEP,
        columns = MANDATORY_COLUMNS,
        new_columns=NEW_COLUMN_NAMES,
        log10p_threshold=LOG10P_THRESHOLD,
        verbose=False,
    ):

    # Preparing result filename keeping only significant QTLs
    res_csv_fname = f"{chr_gz_fname.split('.')[0].split('/')[-1]}_significant_qtls.csv"

    # Don't run again if result file already exists
    if os.path.isfile(f"{res_csv_fname}"):

        log_chr = {"skipped": True}

        if verbose:
            print(
                f"Results for chromosome file {chr_gz_fname} already created at {res_csv_fname}. Skipping.",
                flush=True
            )

    else:

        if verbose:
            print(f"Processing chromosome file:     {chr_gz_fname}.")

        t_start = time.time()

        # Temporarily extracting the .gz file to read its content
        chr_file_gz = protein_tar_file.extractfile(chr_gz_fname)

        # Read .regenie_file and keep only significant QTLs
        summary_stats_significant, log_chr = keep_significant_qtls_from_chr_gz_file(
            chr_file_gz=chr_file_gz,
            separator=separator,
            columns=columns,
            new_columns=new_columns,
            log10p_threshold=log10p_threshold,
            verbose=verbose,
        )

        log_chr["skipped"] = False

        # Saving the significant summary statistics to a new file
        summary_stats_significant.write_csv(f"{res_csv_fname}")
        t_end = time.time()
        if verbose:
            print(f"Significant QTLs saved to: {res_csv_fname}", flush=True)
            print(f"Processed chromosome file {chr_gz_fname} in {t_end - t_start:.2f} s", flush=True)

    return res_csv_fname, log_chr


def merge_significant_qtls_from_all_chr_files(
        csv_fnames,
        output_fname,
        verbose=False,
    ):

    # Read all the significant QTLs from the different chromosome files and concatenate them into one dataframe
    all_significant_qtls = pl.concat(
        [pl.read_csv(csv_fname) for csv_fname in csv_fnames]
    )

    n_kept_qtls = pl.len(all_significant_qtls)

    # Save the concatenated dataframe to a new file
    all_significant_qtls.write_csv(output_fname)
    if verbose:
        print(f"Merged significant QTLs from {len(csv_fnames)} files")
        print(f"All significant QTLs saved to: {output_fname}", flush=True)
        print(f"Total significant QTLs: {n_kept_qtls}", flush=True)

    return all_significant_qtls

def process_one_tar_file(
        tar_fname,
        separator=REGENIE_SEP,
        columns = MANDATORY_COLUMNS,
        new_columns=NEW_COLUMN_NAMES,
        log10p_threshold=LOG10P_THRESHOLD,
        verbose=False,
    ):
    if verbose:
        print(f"Processing protein file: {tar_fname}")

    # There is one .tar file per protein
    with tarfile.open(tar_fname) as protein_tf:

        # In each .tar file, there is one .gz file per chromosome
        list_of_chr_files = protein_tf.getnames()

        log_tar = {
            "tot_chr_files": len(list_of_chr_files),
            "skipped_chr_files": [],
            "n_processed_qtls": 0,
            "processed_qtls": [],
            "log10p_threshold": log10p_threshold,
        }

        # Keep only actual gz files, not the directories
        list_of_chr_files = [f for f in list_of_chr_files if f.endswith('.gz')]

        # Focusing on one chromosome at a time
        # For each chr, there one .gz file containing one .regenie file
        all_csv_fnames = []

        for chr_gz_fname in list_of_chr_files:

            res_csv_fname, log_chr = process_one_chr_from_protein_tar_file(
                protein_tar_file=protein_tf,
                chr_gz_fname=chr_gz_fname,
                separator=separator,
                columns=columns,
                new_columns=new_columns,
                log10p_threshold=log10p_threshold,
                verbose=verbose,
            )

            if log_chr["skipped"]:
                log_tar["skipped_chr_files"].append(chr_gz_fname)
            else:
                log_tar["n_processed_qtls"] += log_chr["n_tot_qtls"]
                log_tar["processed_qtls"] += log_chr["all_qtls"]
            all_csv_fnames.append(res_csv_fname)

    return all_csv_fnames, log_tar

def process_one_region_folder(
        synapse_folder_id=REGION_PQTL_DIR,
        download_location=DOWNLOAD_LOCATION,
        login_kwargs=LOGIN_KWARGS,
        regenie_sep = REGENIE_SEP,
        regenie_columns = MANDATORY_COLUMNS,
        csv_columns = NEW_COLUMN_NAMES,
        log10p_threshold = LOG10P_THRESHOLD,
        create_log = False,
        verbose=False,
):

    t_start_region = time.time()

    # ----------- List of protein tar files in Region folder -----------
    # Get a list of (synapse_id, tar_name) for all the protein tar files
    tar_entities = list_available_protein_tar_files(
        synapse_folder_id=synapse_folder_id,
        login_kwargs=login_kwargs,
    )

    tot_tar_files = len(tar_entities)
    tar_skipped = []
    tar_downloaled_skip = []
    tar_processed = []
    all_merged_csv_fnames = []

    if verbose:
        print(f"Found {tot_tar_files} protein tar files in Synapse folder {synapse_folder_id}.")


    # --------------  Process each protein one by one -----------------
    for synapse_id, tar_name in tar_entities:

        if verbose:
            print(f"Processing Synapse ID: {synapse_id}")

        # Preparing result filename keeping only significant QTLs
        res_merged_fname = f"{tar_name.split('_')[0].split('/')[-1]}_significant_qtls"


        # -------- Skipping if results pre-exist --------------

        # Don't run again if result file already exists
        if os.path.isfile(f"{res_merged_fname}.csv"):
            tar_skipped.append(tar_name)
            if verbose:
                print(
                    f"Result file already created: {res_merged_fname}.csv. Skipping.",
                    flush=True
                )

        else:

            t_start = time.time()

            # ------------ Download tar file --------------

            tar_fname, skipped = download_protein_tar_file(
                synapse_id= synapse_id,
                download_location=download_location,
                login_kwargs=login_kwargs,
                verbose=verbose,
            )

            if skipped:
                tar_downloaled_skip.append(tar_name)

            # ------------ Process tar file --------------

            all_csv_fnames, log_tar = process_one_tar_file(
                tar_fname,
                separator=regenie_sep,
                columns = regenie_columns,
                new_columns=csv_columns,
                log10p_threshold=log10p_threshold,
                verbose=verbose,
            )

            # ------- Merge csv files that kept only significant QTLs --
            all_significant_qtls = merge_significant_qtls_from_all_chr_files(
                csv_fnames=all_csv_fnames,
                output_fname=f"{res_merged_fname}.csv",
                verbose=verbose,
            )

            n_kept_qtls = pl.len(all_significant_qtls)
            log_tar["n_kept_qtls"] = n_kept_qtls

            t_end = time.time()


            # --------- Create log file for the protein ---------------
            if create_log:
                log_fname = f"{res_merged_fname}-log.json"
                log_dict = log_tar + {
                    "synapse_id": synapse_id,
                    "tar_name": tar_name,
                    "all_csv_fnames": all_csv_fnames,
                    "res_merged_csv_fname": f"{res_merged_fname}.csv",
                    "time_taken" : t_end - t_start,
                }
                create_log(log_fname, log_dict, overwrite=True, verbose=verbose)

            tar_processed.append(tar_fname)
            all_merged_csv_fnames.append(f"{res_merged_fname}.csv")

            # --------------- Warnings and verbose ------------------

            if log_tar["skipped_chr_files"]:
                print(f"[WARNING]: Skipped {len(log_tar['skipped_chr_files'])} chromosome files in tar file {tar_name}!")
                print(CONSISTENCY_MESSAGE)

            if verbose:
                print(f"Processed tar file {tar_fname} in {t_end - t_start:.2f} s", flush=True)


    # --------------- Merge tar merged CSV files ------------------
    # Add one column with the protein name

    # TODO!

    t_end_region = time.time()

    # ---------------Create log file for the region ------------------
    if create_log:
        # TODO!
        pass

    if verbose:
        print(f"Processed region tar files in {t_end_region - t_start_region:.2f} s", flush=True)

    if tar_downloaled_skip and verbose:
        print(f"{len(tar_downloaled_skip)}/{tot_tar_files} tar files were already downloaded.")

    if tar_skipped:
        print(f"[WARNING]: {len(tar_skipped)}/{tot_tar_files}  skipped tar files! ")
        print(CONSISTENCY_MESSAGE)

    return all_merged_csv_fnames


# TODO: One function that detects whether the ID is a chr, a tar or a folder and redirected as appropriate

# TODO: Write kept QTLS and those that are not kept
# Merge also the kept QTLs list files (but don't make it mandatory)

# Add info about n_skipped, n_processed regarding the files
# Have a create_log parameter, set to False by default, that create a log file with the info about n_skipped, n_processed, fnames, the time taken for each step, date of execution. The log file will be a CSV file.

# TODO: tests

# TODO: process map files as well, in another file
# -
