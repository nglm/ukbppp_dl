# UKB PPP Download

[![PyPI version](https://badge.fury.io/py/ukbppp-dl.svg)](https://badge.fury.io/py/ukbppp-dl)
[![Documentation Status](https://readthedocs.org/projects/ukbppp-dl/badge/?version=latest)](https://ukbppp-dl.readthedocs.io/en/latest/?badge=latest)


`ukbppp_dl` is a small Python package for easy, robust, memory-efficient and traceable downloading of [UK Biobank Pharma Proteomics Project (UKB-PPP)](https://www.synapse.org/Synapse:syn51364943/wiki/622119) files from Synapse, including pGWAS/pQTL summary statistics with potential filtering based on your given significance threshold.

The most important function is `keep_significant_qtls_from_region`, which can:

- Find all protein tar archives in a given ancestry group in the UKB PPP pGWAS Synapse folder.
- Download only the protein tar files that are needed, and can delete them on the fly to save space.
- Extract chromosome-level REGENIE results from each tar archive.
- Filter rows by a `LOG10P` significance threshold.
- Can write per-chromosome, per-protein, and region-level outputs. Outputs include:
  - CSV files with the significant QTLs only.
  - JSON logs with the run parameters, kept proteins, and output file paths to ensure traceability.
  - A text output file that captures the run transcript.
- Automatically detect and reuse compatible partial results from interrupted runs (because you needed to restart your local machine or your remote server shut down, etc.).
- Clean up intermediate files when requested.

See section [Basic usage](https://ukbppp-dl.readthedocs.io/en/latest/index.html#basic-usage) for more practical details.

## Why this package matters

Main strengths of the package:

1. **Memory friendly**: Instead of downloading the whole region file (probably around 9PT), `ukbppp_dl` processes protein individually, and can delete intermediate files on the fly, making it memory-friendly. In addition it doesn't require to extract the whole archives to the disk and read the chromosome files directly from the tar archives instead. You won't need much more than 1GB of free space to run the function.
2. **Restart-friendly**: compatible partial results, logs and outputs can be reused instead of recomputed. Interrupt your run without worrying about losing all the progress, nor creating inconsistent results.
3. **Traceability and reproducibility**: It keeps the results traceable through JSON logs and text output files. Thus, once the final CSV files are generated, you can still check with which parameters they were generated.
4. **Flexible cleanup control**: It gives you control over cleanup, generated files and verbosity, so that the final outputs are easy to understand and manage, avoiding much manual intervention (and thus potential errors and inconsistencies).
5. **Practical**: No need to go into the details of Synapse client, tar file handling, and data processing. The main function can be used with a single call, and the parameters are intuitive.

## Documentation

The full documentation is available at [UKB PPP Download - Documentation](https://ukbppp-dl.readthedocs.io/).

## Installation

### Alternative 1: With uv:

```bash
# From PyPI
uv add ukbppp-dl
# Alternatively, from github directly
uv add "ukbppp-dl @ git+https://github.com/nglm/ukbppp-dl.git"
```

### Alternative 2: With poetry:

```bash
# From PyPI
poetry add ukbppp-dl
# Alternatively, from github directly
poetry add git+https://github.com/nglm/ukbppp-dl.git
```

### Alternative 3: With pip:

```bash
# From PyPI
pip install ukbppp-dl
# Alternatively, from github directly
pip install git+https://github.com/nglm/ukbppp_dl.git
```

## Requirements

You need a Synapse account with access to [UK Biobank Pharma Proteomics Project (UKB-PPP)](https://www.synapse.org/Synapse:syn51364943/wiki/622119). `ukbppp_dl` can use your local Synapse configuration file `~/.synapseConfig` if you already have one set up. Alternatively, you can also provide your Synapse credentials when calling the functions.

An example of the content of the Synapse configuration file `~/.synapseConfig` is shown below. Replace the placeholders with your actual Synapse credentials. You can find more details about the Synapse configuration file in the [Synapse documentation](https://python-docs.synapse.org/en/stable/tutorials/authentication/).

```yaml
[default]
username = your.synapse.email@mail.org
authtoken = YouRAuthenticationKeyWithManyCharacters

[cache]
location = ~/.synapseCache
```

## Basic usage

You can see an example script in `scripts/download_significant_qtls.py` and adapt it to your needs. The main function to use is `keep_significant_qtls_from_region`, which can be used as follows:

```python
from ukbppp_dl.pgwas import keep_significant_qtls_from_region, PGWAS_REGIONS

# Synapse directory containing pQTL summary statistics (here for Europe)
REGION = PGWAS_REGIONS["European"]

# Significance threshold for pQTLs (LOG10P > 7 corresponds to p-value < 1e-7)
LOG10P_THRESHOLD = 7

# Whether to create a log file
# (0: no log file, >0: create different levels of log files)
CREATE_LOG = 2

# Whether to have an output text describing the function's run
# (0: no text, >0: create different levels of verbosity)
VERBOSE = 3

# set to a list of protein tar file names or synapse IDs if you want to process only specific proteins
# PROTEIN_TO_PROCESS = ["ACOT13_Q9NPJ3_OID31522_v1_Oncology_II.tar", "syn52363271"]

# otherwise set to None to process all proteins in the region
PROTEIN_TO_PROCESS = None

all_significant_qtls, log_reg = keep_significant_qtls_from_region(
	synapse_folder_id=REGION,
	download_location="./data",
	res_location="./results",
	log10p_threshold=LOG10P_THRESHOLD,
	create_log=CREATE_LOG,
	verbose=VERBOSE,
	delete_downloaded_tar=True,
	delete_chr_csv=True,
    protein_to_process=PROTEIN_TO_PROCESS,
	delete_tar_csv=False,
	delete_tar_log=False,
	delete_partial_logs=False,
	delete_partial_outputs=False,
)
```

The resulting `all_significant_qtls` is a [Polars DataFrame](https://docs.pola.rs/py-polars/html/reference/dataframe/index.html) containing all significant QTLs across the selected region. `log_reg` is the final region log dictionary, which records the parameters, processed protein tar files, kept proteins, and output file paths, etc. used when producing `all_significant_qtls`.

See an example of resulting `all_significant_qtls` dataframe below (with a subset of columns, and fictive BETA, SE and LOG10P values):

| Protein_name | CHROM | POS       | ID                     | BETA  | SE   | LOG10P |
| ------------ | ----- | --------- | ---------------------- | ----- | ---- | ------ |
| ABCA2        | 11    | 1737296   | 11:1758526:C:T:imp:v1  | -0.5  | 0.1  | 8.0    |
| ABCA2        | 11    | 1764171   | 11:1785401:T:C:imp:v1  | -0.4  | 0.2  | 7.5    |
| ABCA2        | 9     | 137017531 | 9:139911983:T:G:imp:v1 | -0.09 | 0.01 | 13.6   |
| ABCA2        | 9     | 137017726 | 9:139912178:G:A:imp:v1 | 0.05  | 0.03 | 8.1    |
| ABHD14B      | 10    | 63118358  | 10:64878118:C:G:imp:v1 | 0.05  | 0.01 | 8.5    |
| ABHD14B      | 10    | 63122540  | 10:64882300:C:G:imp:v1 | 0.04  | 0.02 | 8.4    |
| ABHD14B      | 12    | 54342686  | 12:54736470:A:G:imp:v1 | 0.04  | 0.03 | 8.1    |
| ABHD14B      | 18    | 57633880  | 18:55301112:A:G:imp:v1 | -0.06 | 0.01 | 7.3    |
| ABHD14B      | 19    | 55014977  | 19:55526345:T:G:imp:v1 | 0.06  | 0.01 | 7.2    |
| ABHD14B      | 19    | 55025227  | 19:55536595:G:A:imp:v1 | 0.07  | 0.01 | 8.0    |
| ...          | ...   | ...       | ...                    | ...   | ...  | ...    |


## Contribute

- Issue Tracker: [github.com/nglm/ukbppp_dl/issues](https://github.com/nglm/ukbppp_dl/issues).
- Source Code: [github.com/nglm/ukbppp_dl](github.com/nglm/ukbppp_dl).

## Support

If you are having issues, [please let me know](https://www.uib.no/en/persons/Natacha.Madeleine.Georgette.Galmiche) or [create an issue](https://github.com/nglm/ukbppp_dl/issues).

## License

MIT License. See `LICENSE.txt`.
