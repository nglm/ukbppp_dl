import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)

import tarfile
import gzip
import polars as pl
import os
import json
import pathlib


from ..pqtls import (
       list_available_protein_tar_files, download_protein_tar_file,
       process_one_chr_from_protein_tar_file,
       process_one_tar_file, merge_significant_qtls_from_all_chr_files,
       process_one_region_folder,
       save_log,
    )

# Synapse directory containing pQTL summary statistics (here for Europe)
REGION_PQTL_DIR = 'syn51365303'

DOWNLOAD_LOCATION = "./data"
RES_LOCATION = "./results_tests"

# Specific protein tar files for test purpose
ACOT13_ID = "syn52362654"
ACOT13_FNAME = "ACOT13_Q9NPJ3_OID31522_v1_Oncology_II.tar"
ZNF174_ID = "syn52363271"

additional_protein_ids = [
    "syn52361344",
    "syn51470065",
    "syn52363381",
    "syn51470344",
    "syn52363597",
]

# Path to files
ACOT13_PATH = f"{DOWNLOAD_LOCATION}/{ACOT13_FNAME}"

# Specific chr files for test purpose
ACOT13_CHR1_FILE = "discovery_chr1_ACOT13:Q9NPJ3:OID31522:v1:Oncology_II.gz"

# Mandatory columns in the .regenie files
MANDATORY_COLUMNS = ["CHROM", "GENPOS", "ID", "BETA", "SE", "LOG10P"]
NEW_COLUMN_NAMES = ["chrom", "qtl_pos", "qtl_id", "beta", "se", "log10p"]

# Significance threshold for pQTLs
LOG10P_THRESHOLD = 7

# Separator used in the .regenie files
REGENIE_SEP = " "

# Whether to create a log file
CREATE_LOG = False

LOG_KWARGS = {
    "overwrite": False,
    "add_date": True,
    "new_name": False,
    "verbose": True,
}
# Synapse login kwargs
LOGIN_KWARGS = {}

# Actual number of tar files in the directory
N_EXPECTED_TAR_FILES = 2940

# Actual number of QTLs in chr 1 file for ACOT13
N_EXPECTED_QTLS_ACOT13_CHR1 = 1212345

# Actual number of QTLs in ACOT13
N_EXPECTED_QTLS_ACOT13 = 16025237

# Actual number of chromosome files in the tar file
N_EXPECTED_CHR_FILES = 23

def test_save_log():
    dict = {"key1": "value", "key2": 123, "key3": [1, 2, 3]}
    output_name = f"{RES_LOCATION}/test-log.json"

    # Creating new log file
    log_fname = save_log(
        output_name,
        dict,
        overwrite=True,
        add_date=False,
    )

    assert isinstance(log_fname, str)
    assert log_fname == output_name
    assert os.path.isfile(log_fname)
    with open(log_fname, "r") as f:
        log_content = json.load(f)
    assert log_content == dict

    # Trying to create log file but actually skip
    log_fname = save_log(
        output_name,
        dict,
        overwrite=False,
        new_name=False,
        add_date=False,
    )
    assert log_fname is None
    assert os.path.isfile(output_name)
    with open(output_name, "r") as f:
        log_content = json.load(f)
    assert log_content == dict

    # Trying to create log file but use new name
    log_fname = save_log(
        output_name,
        dict,
        overwrite=False,
        new_name=True,
        add_date=False,
    )
    assert isinstance(log_fname, str)
    assert log_fname != output_name
    assert os.path.isfile(output_name)
    assert os.path.isfile(log_fname)
    with open(log_fname, "r") as f:
        log_content = json.load(f)
    assert log_content == dict

def test_list_available_protein_tar_files():

        tar_files = list_available_protein_tar_files()

        # Is a list of tuples (id, name)
        assert isinstance(tar_files, list)
        assert len(tar_files) == N_EXPECTED_TAR_FILES
        assert all(len(tuple) == 2 for tuple in tar_files)

        names = [name for _, name in tar_files]
        assert all(isinstance(name, str) for name in names)
        assert all(name.endswith('.tar') for name in names)
        IDs = [id for id, _ in tar_files]
        assert all(isinstance(id, str) for id in IDs)
        assert all(id.startswith('syn') for id in IDs)


def test_download_protein_tar_file():
        warnings.filterwarnings("ignore", category=DeprecationWarning)

        expected_fname, skipped = download_protein_tar_file(
                ACOT13_ID,
                download_location=DOWNLOAD_LOCATION,
                verbose = True,
        )

        assert isinstance(expected_fname, str)
        assert expected_fname.endswith(".tar")
        assert expected_fname == ACOT13_PATH

        assert isinstance(skipped, bool)

def test_process_one_chr_from_protein_tar_file():

    warnings.filterwarnings("ignore", category=DeprecationWarning)

    expected_fname, skipped = download_protein_tar_file(
        ACOT13_ID,
        download_location=DOWNLOAD_LOCATION,
        verbose = True,
    )

    with tarfile.open(ACOT13_PATH) as protein_tf:

        res_csv_fname, log = process_one_chr_from_protein_tar_file(
            protein_tf,
            ACOT13_CHR1_FILE,
            res_location=RES_LOCATION,
            separator=REGENIE_SEP,
            columns = MANDATORY_COLUMNS,
            new_columns=NEW_COLUMN_NAMES,
            log10p_threshold=2,
            create_log=2,
            verbose=True,
        )

    assert isinstance(res_csv_fname, str)
    assert res_csv_fname.endswith(".csv")
    assert os.path.isfile(res_csv_fname)

    assert isinstance(log, dict)
    if log["skipped"]:
        assert True
    else:
        keys = ["log10p_threshold", "n_tot_qtls", "n_kept_qtls", "source_chr_file", "skipped", "log_filename"]
        assert all(key in log for key in keys), f"Missing keys in log: {[key for key in keys if key not in log]}"
        assert log["log10p_threshold"] == 2
        assert log["n_tot_qtls"] >= log["n_kept_qtls"]
        assert log["n_tot_qtls"] == N_EXPECTED_QTLS_ACOT13_CHR1


def test_process_one_tar_file():

    warnings.filterwarnings("ignore", category=DeprecationWarning)

    expected_fname, skipped = download_protein_tar_file(
        ACOT13_ID,
        download_location=DOWNLOAD_LOCATION,
        verbose = True,
    )

    all_csv_fnames, log = process_one_tar_file(
            expected_fname,
            res_location=RES_LOCATION,
            separator=REGENIE_SEP,
            columns = MANDATORY_COLUMNS,
            new_columns=NEW_COLUMN_NAMES,
            log10p_threshold=LOG10P_THRESHOLD,
            create_log=2,
            verbose=True,
    )

    assert isinstance(all_csv_fnames, list)
    assert all(isinstance(fname, str) for fname in all_csv_fnames)
    assert all(fname.endswith(".csv") for fname in all_csv_fnames)
    assert all(os.path.isfile(fname) for fname in all_csv_fnames)

    assert isinstance(log, dict)
    keys = [
        "log_filename", "tar_fname", "protein_name",
        "tot_chr_files", "skipped_chr_files", "all_csv_fnames",
        "n_processed_qtls", "log10p_threshold", "regenie_columns",
        "csv_columns"]
    assert all(key in log for key in keys), f"Missing keys in log: {[key for key in keys if key not in log]}"
    assert log["tot_chr_files"] == N_EXPECTED_CHR_FILES
    assert len(log["skipped_chr_files"]) >= 0
    assert len(log["skipped_chr_files"]) <= N_EXPECTED_CHR_FILES


    if log["skipped_chr_files"] == N_EXPECTED_CHR_FILES:
        assert log["n_processed_qtls"] == 0
    elif log["skipped_chr_files"] == 0:
        assert log["n_processed_qtls"] == N_EXPECTED_QTLS_ACOT13

def test_merge_significant_qtls_from_all_chr_files():
    # This test relies on the previous one, which processes the tar file and creates the csv files
    # If the previous test fails, this one will also fail

    all_csv_fnames, log_tar = process_one_tar_file(
            ACOT13_PATH,
            res_location=RES_LOCATION,
            separator=REGENIE_SEP,
            columns = MANDATORY_COLUMNS,
            new_columns=NEW_COLUMN_NAMES,
            log10p_threshold=LOG10P_THRESHOLD,
            create_log=0,
            verbose=True,
    )

    out_fname = f"{RES_LOCATION}/ACOT13-test_merging_significant_qtls.csv"

    all_significant_qtls, log = merge_significant_qtls_from_all_chr_files(
        all_csv_fnames,
        output_fname=out_fname,
        create_log=2,
        delete_chr_csv=False,
        verbose=True,
    )

    assert isinstance(all_significant_qtls, pl.DataFrame)
    assert isinstance(log, dict)

    keys = ["log_filename", "merged_csv_fname", "n_chr_files_merged", "n_kept_qtls", "n_kept_qtls_per_chr_file", "min_log10p"]
    assert all(key in log for key in keys), f"Missing keys in log: {[key for key in keys if key not in log]}"


def test_process_one_region_folder():

    warnings.filterwarnings("ignore", category=DeprecationWarning)

    all_significant_qtls, log_reg = process_one_region_folder(
        synapse_folder_id=REGION_PQTL_DIR,
        download_location=DOWNLOAD_LOCATION,
        res_location=RES_LOCATION,
        login_kwargs=LOGIN_KWARGS,
        regenie_sep=REGENIE_SEP,
        regenie_columns = MANDATORY_COLUMNS,
        csv_columns=NEW_COLUMN_NAMES,
        log10p_threshold=LOG10P_THRESHOLD,
        create_log=3,
        log_kwargs=LOG_KWARGS,
        protein_to_process = [ACOT13_ID, ZNF174_ID],
        verbose=1,
        delete_downloaded_tar = False,
        delete_chr_csv = True,
        delete_tar_csv = True,
        delete_partial_logs = True,
    )

    assert isinstance(all_significant_qtls, pl.DataFrame)
    assert isinstance(log_reg, dict)


    # In this test we will reuse a part-file (that we hopfully won't delete each time)
    all_significant_qtls, log_reg = process_one_region_folder(
        synapse_folder_id=REGION_PQTL_DIR,
        download_location=DOWNLOAD_LOCATION,
        res_location=RES_LOCATION,
        login_kwargs=LOGIN_KWARGS,
        regenie_sep=REGENIE_SEP,
        regenie_columns = MANDATORY_COLUMNS,
        csv_columns=NEW_COLUMN_NAMES,
        log10p_threshold=LOG10P_THRESHOLD,
        create_log=2,
        log_kwargs=LOG_KWARGS,
        protein_to_process = additional_protein_ids,
        verbose=3,
        delete_downloaded_tar = False,
        delete_chr_csv = True,
        delete_tar_csv = True,
        delete_tar_log = True,
        delete_partial_logs = "current",
    )

    assert isinstance(all_significant_qtls, pl.DataFrame)
    assert isinstance(log_reg, dict)

    all_significant_qtls, log_reg = process_one_region_folder(
        synapse_folder_id=REGION_PQTL_DIR,
        download_location=DOWNLOAD_LOCATION,
        res_location=RES_LOCATION,
        login_kwargs=LOGIN_KWARGS,
        regenie_sep=REGENIE_SEP,
        regenie_columns = MANDATORY_COLUMNS,
        csv_columns=NEW_COLUMN_NAMES,
        log10p_threshold=LOG10P_THRESHOLD,
        create_log=0,
        log_kwargs=LOG_KWARGS,
        protein_to_process = [ACOT13_ID, ZNF174_ID],
        verbose=2,
        delete_downloaded_tar = False,
        delete_chr_csv = True,
        delete_tar_csv = True,
        delete_tar_log = True,
        delete_partial_logs = "all",
    )

    assert isinstance(all_significant_qtls, pl.DataFrame)
    assert isinstance(log_reg, dict)
