
from ..pqtls import process_one_tar_file

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

def test_process_one_tar_file():

        all_csv_fnames, log_tar = process_one_tar_file(
                ACOT13_FILE,
                separator=REGENIE_SEP,
                columns = MANDATORY_COLUMNS,
                new_columns=NEW_COLUMN_NAMES,
                log10p_threshold=LOG10P_THRESHOLD,
                verbose=False,
        )
