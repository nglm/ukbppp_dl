from ukbppp_dl.pgwas import keep_significant_qtls_from_region, PGWAS_REGIONS


# Synapse directory containing pQTL summary statistics (here for Combined)
REGION = PGWAS_REGIONS["Combined"]

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
# (0: no log file, >0: create different levels of log files)
CREATE_LOG = 2

# Whether to have an output text describing the function's run
# (0: no text, >0: create different levels of verbosity)
VERBOSE = 3


all_significant_qtls, log_reg = keep_significant_qtls_from_region(
        synapse_folder_id=REGION,
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