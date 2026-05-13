import warnings
import os
import json

from ..common import ( download_from_synapse, save_log, )

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
NEW_COLUMN_NAMES = ["CHR", "QTL_POS", "QTL_ID", "BETA", "SE", "LOG10P"]

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


def test_download_from_synapse():
        warnings.filterwarnings("ignore", category=DeprecationWarning)

        expected_fname, skipped = download_from_synapse(
                ACOT13_ID,
                download_location=DOWNLOAD_LOCATION,
                verbose = True,
        )

        assert isinstance(expected_fname, str)
        assert expected_fname.endswith(".tar")
        assert expected_fname == ACOT13_PATH

        assert isinstance(skipped, bool)
