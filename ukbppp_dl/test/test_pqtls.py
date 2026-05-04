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
       process_one_tar_file,
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

# Synapse login kwargs
LOGIN_KWARGS = {}

# Actual number of tar files in the directory
N_EXPECTED_TAR_FILES = 2940

# Actual number of QTLs in chr 1 file for ACOT13
N_EXPECTED_QTLS_CHR1 = 1212345

# Actual number of QTLs in chr 1 file for ACOT13
N_EXPECTED_QTLS = 16025237

# Actual number of chromosome files in the tar file
N_EXPECTED_CHR_FILES = 23

def test_save_log():
    is_new = save_log(
        "test-log.json",
        {"key1": "value", "key2": 123, "key3": [1, 2, 3]},
        overwrite=True
    )

    assert is_new == True

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

        res_csv_fname, log_chr = process_one_chr_from_protein_tar_file(
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

    assert isinstance(log_chr, dict)
    if log_chr["skipped"]:
        assert True
    else:
        keys = ["log10p_threshold", "n_tot_qtls", "n_kept_qtls", "source_chr_file", "skipped"]
        assert all(key in log_chr for key in keys)
        assert log_chr["log10p_threshold"] == 2
        assert log_chr["n_tot_qtls"] >= log_chr["n_kept_qtls"]
        assert len(log_chr["all_qtls"]) == log_chr["n_tot_qtls"]
        assert log_chr["n_tot_qtls"] == N_EXPECTED_QTLS_CHR1


def test_process_one_tar_file():

    warnings.filterwarnings("ignore", category=DeprecationWarning)

    expected_fname, skipped = download_protein_tar_file(
        ACOT13_ID,
        download_location=DOWNLOAD_LOCATION,
        verbose = True,
    )

    all_csv_fnames, log_tar = process_one_tar_file(
            ACOT13_PATH,
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

    assert isinstance(log_tar, dict)
    keys = ["tot_chr_files", "skipped_chr_files", "n_processed_qtls", "log10p_threshold"]
    assert all(key in log_tar for key in keys)
    assert log_tar["tot_chr_files"] == N_EXPECTED_CHR_FILES
    assert len(log_tar["skipped_chr_files"]) >= 0
    assert len(log_tar["skipped_chr_files"]) < N_EXPECTED_CHR_FILES


    if log_tar["skipped_chr_files"] == N_EXPECTED_CHR_FILES:
        assert log_tar["n_processed_qtls"] == 0
    elif log_tar["skipped_chr_files"] == 0:
        assert log_tar["n_processed_qtls"] == N_EXPECTED_QTLS

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

    

