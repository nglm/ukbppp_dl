from ukbppp_dl.pgwas import process_one_region_folder


# Synapse directory containing pQTL summary statistics (here for Europe)
REGION_PQTL_DIR = 'syn51365303'

DOWNLOAD_LOCATION = "./data"
RES_LOCATION = "./results"

# Mandatory columns in the .regenie files
MANDATORY_COLUMNS = ["CHROM", "GENPOS", "ID", "BETA", "SE", "LOG10P"]
# New column names to use in the output csv files
# Note that in any case "LOG10P" must be kept and with the same name
NEW_COLUMN_NAMES = ["CHR", "QTL_POS", "QTL_ID", "BETA", "SE", "LOG10P"]

# Significance threshold for pQTLs
LOG10P_THRESHOLD = 7

# Whether to create a log file
CREATE_LOG = 2

VERBOSE = 3


all_significant_qtls, log_reg = process_one_region_folder(
        synapse_folder_id=REGION_PQTL_DIR,
        download_location=DOWNLOAD_LOCATION,
        res_location=RES_LOCATION,
        regenie_columns = MANDATORY_COLUMNS,
        csv_columns=NEW_COLUMN_NAMES,
        log10p_threshold=LOG10P_THRESHOLD,
        create_log=CREATE_LOG,
        verbose=VERBOSE,
        delete_downloaded_tar = True,
        delete_chr_csv = True,
        delete_tar_csv = False,
        delete_tar_log = False,
        delete_partial_logs = False,
        delete_partial_outputs = False,
    )