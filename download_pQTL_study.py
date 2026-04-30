import time

import synapseclient
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

# Synapse login kwargs
LOGIN_KWARGS = {}

def download_protein_tar_file(
        synapse_id,
        download_location=".",
        login_kwargs={}
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
    print(f"Expected file path: {expected_file_path}")
    if os.path.isfile(expected_file_path):
        print(f"File already downloaded: {expected_file_path}. Skipping download.")
    else:
        # Download the file to the specified location
        entity = syn.get(synapse_id, downloadLocation=download_location, downloadFile=True)
        expected_file_path = entity.path
        print(f"Downloaded file: {expected_file_path}")

    return expected_file_path

def keep_significant_qtls_from_chr_gz_file(
        chr_file_gz,
        separator=REGENIE_SEP,
        columns = MANDATORY_COLUMNS,
        new_columns=NEW_COLUMN_NAMES,
        log10p_threshold=LOG10P_THRESHOLD,
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
    ):
    print(f"Processing chromosome file: {chr_gz_fname}")

    # Preparing result filename keeping only significant QTLs
    res_csv_fname = f"{chr_gz_fname.split('.')[0].split('/')[-1]}_significant_qtls.csv"

    # Don't run again if result file already exists
    if os.path.isfile(f"{res_csv_fname}"):
        print(
            f"Result file already created: {res_csv_fname}. Skipping.",
            flush=True
        )
        pass

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
        )

        # Saving the significant summary statistics to a new file
        summary_stats_significant.write_csv(f"{res_csv_fname}")
        print(f"Significant QTLs saved to: {res_csv_fname}", flush=True)
        t_end = time.time()
        print(f"Processed chromosome file: {chr_gz_fname} in {t_end - t_start:.2f} s", flush=True)

    return res_csv_fname


def process_one_tar_file(
        tar_full_fname,
        separator=REGENIE_SEP,
        columns = MANDATORY_COLUMNS,
        new_columns=NEW_COLUMN_NAMES,
        log10p_threshold=LOG10P_THRESHOLD,
    ):

    # There is one .tar file per protein
    with tarfile.open(tar_full_fname) as protein_tf:

        # In each .tar file, there is one .gz file per chromosome
        list_of_chr_files = protein_tf.getnames()


        # print('Code executed in %.2f s' %(t_end - t_start))
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
            )

            all_csv_fnames.append(res_csv_fname)

        t_end = time.time()
        print(f"Processed protein file {tar_full_fname} in {t_end - t_start:.2f} s", flush=True)

        return all_csv_fnames


def merge_significant_qtls_from_all_chr_files(
        csv_fnames,
        output_fname,
    ):

    # Read all the significant QTLs from the different chromosome files and concatenate them into one dataframe
    all_significant_qtls = pl.concat(
        [pl.read_csv(csv_fname) for csv_fname in csv_fnames]
    )

    # Save the concatenated dataframe to a new file
    all_significant_qtls.write_csv(output_fname)
    print(f"All significant QTLs saved to: {output_fname}", flush=True)

    return output_fname

tar_full_fname = download_protein_tar_file(
    ZNF174_FILE,
    download_location=DOWNLOAD_LOCATION, login_kwargs=LOGIN_KWARGS
)


# TODO: function to list available tar files

# TODO: Do something when the full tar file is processed so that we know when to skip a tar file.

