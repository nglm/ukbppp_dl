import time

import synapseclient
from synapseclient.models import Folder
import tarfile
import gzip
import polars as pl
import os


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

def list_available_protein_tar_files(
        synapse_id=REGION_PQTL_DIR,
        login_kwargs={},
        create_log=False,
        verbose=False,
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
        create_log=False,
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
    if verbose:
        print(f"Expected file path: {expected_file_path}")
    if os.path.isfile(expected_file_path):
        if verbose:
            print(f"File already downloaded: {expected_file_path}. Skipping download.")
    else:
        # Download the file to the specified location
        entity = syn.get(
            synapse_id, downloadLocation=download_location, downloadFile=True
        )
        expected_file_path = entity.path
        if verbose:
            print(f"Downloaded file: {expected_file_path}")

    return expected_file_path

def keep_significant_qtls_from_chr_gz_file(
        chr_file_gz,
        separator=REGENIE_SEP,
        columns = MANDATORY_COLUMNS,
        new_columns=NEW_COLUMN_NAMES,
        log10p_threshold=LOG10P_THRESHOLD,
        create_log=False,
        verbose=False,
    ):
    """Read a .regenie file from a .gz file and keep only significant QTLs.

    Parameters:
    chr_file_gz (file-like object): The .gz file containing the .regenie file for a chromosome.
    separator (str): The separator used in the .regenie file.
    columns (list): The list of mandatory columns to read from the .regenie file.
    new_columns (list): The list of new column names to assign to the dataframe.
    log10p_threshold (float): The threshold for log10p to consider a QTL as significant.

    Returns:
    pl.DataFrame: A dataframe containing only the significant QTLs.
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

        # Filtering the summary statistics to keep only significant QTLs
        summary_stats_significant = summary_stats.filter(
            pl.col("log10p") > log10p_threshold
        )

    return summary_stats_significant

def process_one_chr_from_protein_tar_file(
        protein_tar_file,
        chr_gz_fname,
        separator=REGENIE_SEP,
        columns = MANDATORY_COLUMNS,
        new_columns=NEW_COLUMN_NAMES,
        log10p_threshold=LOG10P_THRESHOLD,
        create_log=False,
        verbose=False,
    ):
    if verbose:
        print(f"Processing chromosome file: {chr_gz_fname}")

    # Preparing result filename keeping only significant QTLs
    res_csv_fname = f"{chr_gz_fname.split('.')[0].split('/')[-1]}_significant_qtls.csv"

    # Don't run again if result file already exists
    if os.path.isfile(f"{res_csv_fname}"):

        if verbose:
            print(
                f"Result file already created: {res_csv_fname}. Skipping.",
                flush=True
            )

    else:

        t_start = time.time()

        # Temporarily extracting the .gz file to read its content
        chr_file_gz = protein_tar_file.extractfile(chr_gz_fname)

        # Read .regenie_file and keep only significant QTLs
        summary_stats_significant = keep_significant_qtls_from_chr_gz_file(
            chr_file_gz=chr_file_gz,
            separator=separator,
            columns=columns,
            new_columns=new_columns,
            log10p_threshold=log10p_threshold,
            create_log=create_log,
            verbose=verbose,
        )

        # Saving the significant summary statistics to a new file
        summary_stats_significant.write_csv(f"{res_csv_fname}")
        t_end = time.time()
        if verbose:
            print(f"Significant QTLs saved to: {res_csv_fname}", flush=True)
            print(f"Processed chromosome file: {chr_gz_fname} in {t_end - t_start:.2f} s", flush=True)

    return res_csv_fname


def merge_significant_qtls_from_all_chr_files(
        csv_fnames,
        output_fname,
        create_log=False,
        verbose=False,
    ):

    # Read all the significant QTLs from the different chromosome files and concatenate them into one dataframe
    all_significant_qtls = pl.concat(
        [pl.read_csv(csv_fname) for csv_fname in csv_fnames]
    )

    # Save the concatenated dataframe to a new file
    all_significant_qtls.write_csv(output_fname)
    if verbose:
        print(f"All significant QTLs saved to: {output_fname}", flush=True)

    return output_fname

def process_one_tar_file(
        tar_fname,
        separator=REGENIE_SEP,
        columns = MANDATORY_COLUMNS,
        new_columns=NEW_COLUMN_NAMES,
        log10p_threshold=LOG10P_THRESHOLD,
        create_log=False,
        verbose=False,
    ):
    if verbose:
        print(f"Processing protein file: {tar_fname}")

    # There is one .tar file per protein
    with tarfile.open(tar_fname) as protein_tf:

        # In each .tar file, there is one .gz file per chromosome
        list_of_chr_files = protein_tf.getnames()

        if verbose:
            print(f"One .gz file per chromosome in the tar file: {list_of_chr_files}")
        t_start = time.time()

        # Keep only actual gz files, not the directories
        list_of_chr_files = [f for f in list_of_chr_files if f.endswith('.gz')]

        # Focusing on one chromosome at a time
        # For each chr, there one .gz file containing one .regenie file
        all_csv_fnames = []
        for chr_gz_fname in list_of_chr_files:

            res_csv_fname = process_one_chr_from_protein_tar_file(
                protein_tar_file=protein_tf,
                chr_gz_fname=chr_gz_fname,
                separator=separator,
                columns=columns,
                new_columns=new_columns,
                log10p_threshold=log10p_threshold,
                create_log=create_log,
                verbose=verbose,
            )

            all_csv_fnames.append(res_csv_fname)

        t_end = time.time()
        if verbose:
            print(f"Processed protein file {tar_fname} in {t_end - t_start:.2f} s", flush=True)

    return all_csv_fnames

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

    # Get a list of (synapse_id, tar_name) for all the protein tar files
    tar_entities = list_available_protein_tar_files(
        synapse_folder_id=synapse_folder_id,
        login_kwargs=login_kwargs,
        create_log=create_log,
        verbose=verbose,
    )

    all_merged_csv_fnames = []

    # Process each protein one by one
    for synapse_id, tar_name in tar_entities:

        if verbose:
            print(f"Processing Synapse ID: {synapse_id}")

        # Preparing result filename keeping only significant QTLs
        res_merged_csv_fname = f"{tar_name.split('_')[0][-1]}_significant_qtls.csv"

        # Don't run again if result file already exists
        if os.path.isfile(f"{res_merged_csv_fname}"):
            if verbose:
                print(
                    f"Result file already created: {res_merged_csv_fname}. Skipping.",
                    flush=True
                )

        else:

            tar_fname = download_protein_tar_file(
                synapse_id= synapse_id,
                download_location=download_location,
                login_kwargs=login_kwargs,
                create_log=create_log,
                verbose=verbose,
            )


            all_csv_fnames = process_one_tar_file(
                tar_fname,
                separator=regenie_sep,
                columns = regenie_columns,
                new_columns=csv_columns,
                log10p_threshold=log10p_threshold,
                create_log=create_log,
                verbose=verbose,
            )

            # Now merge all csv files that kept only significant QTLs
            res_merged_csv_fname = merge_significant_qtls_from_all_chr_files(
                csv_fnames=all_csv_fnames,
                output_fname=res_merged_csv_fname,
                create_log=create_log,
                verbose=verbose,
            )

            all_merged_csv_fnames.append(res_merged_csv_fname)

    return all_merged_csv_fnames


# TODO: One function that detects whether the ID is a chr, a tar or a folder and redirected as appropriate

# TODO: Write kept QTLS and those that are not kept
# Merge also the kept QTLs list files (but don't make it mandatory)

# Add info about n_skipped, n_processed regarding the files
# Have a create_log parameter, set to False by default, that create a log file with the info about n_skipped, n_processed, fnames, the time taken for each step, date of execution. The log file will be a CSV file.

# TODO: tests

# TODO: process map files as well, in another file
# -
