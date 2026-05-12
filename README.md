# UKB PPP Download

`ukbppp_dl` is a small Python package for (easy, robust, memory-efficient and traceable) downloading and filtering UK Biobank PPP pGWAS/pQTL summary statistics from Synapse.

The most important function is `ukbppp_dl.pgwas.keep_significant_qtls_from_region`, which can:

- Finds all protein tar archives in a given ancestry group in the UKB PPP pGWAS Synapse folder.
- Download only the protein tar files that are needed, and can delete them on the fly to save space.
- Extract chromosome-level REGENIE results from each tar archive
- Filter rows by a `LOG10P` significance threshold
- Can write per-chromosome, per-protein, and region-level outputs. Outputs include:
  - CSV files with the significant QTLs only
  - JSON logs with the run parameters, kept proteins, and output file paths to ensure traceability
  - a text output file that captures the run transcript
- Automatically detect and reuse compatible partial results from interrupted runs.
- Clean up intermediate files when requested

## Why this package matters

Main strengths of the package:

1. **Memory friendly**: Instead of downloading the whole region file (probably around 9PT), `ukbppp_dl` processes protein individually, and can delete intermediate files on the fly, making it memory-friendly. In addition it doesn't require to extract the whole archives to the disk and read the chromosome files directly from the tar archives instead. You won't need much more than 1GB of free space to run the function.
2. **Restart-friendly**: compatible partial logs and outputs can be reused instead of recomputed.
3. **Traceability and reproducibility**: It keeps the results traceable through JSON logs and text output files. Once the final CSV files are generated, one can still check with which parameters they were generated.
4. **Flexible cleanup control**: It gives you control over cleanup, generated files and verbosity, so that the final outputs are easy to understand and manage, avoiding manual intervention (and potential errors).
5. **Practical**: No need to go into the details of Synapse client, tar file handling, and data processing. The main function can be used with a single call, and the parameters are intuitive.

## Documentation

The full documentation is available at XXX.

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
pip install git+https://github.com/nglm/ukbppp-dl.git
```

### Alternative 4: With anaconda:

```bash
# activate your environment (replace myEnv with your environment name)
conda activate myEnv
# install pip first in your environment
conda install pip
# install ukbppp-dl on your anaconda environment with pip
pip install ukbppp-dl
```

## Requirements

You need a Synapse account with access to the target folder or file IDs. `ukbppp_dl` can use your local Synapse configuration file `~/.synapseConfig` if you already have one set up. Alternatively, you can also provide your Synapse credentials when calling the functions.

An example of the content of the Synapse configuration file `~/.synapseConfig` is shown below. You can find more details about the Synapse configuration file in the [Synapse documentation](https://python-docs.synapse.org/en/stable/tutorials/authentication/).

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
	download_location="./data",
	res_location="./results",
	log10p_threshold=7,
	create_log=CREATE_LOG,
	verbose=VERBOSE,
	delete_downloaded_tar=True,
	delete_chr_csv=True,
	delete_tar_csv=False,
	delete_tar_log=False,
	delete_partial_logs=False,
	delete_partial_outputs=False,
)
```

The resulting `all_significant_qtls` is a Polars DataFrame containing all significant QTLs across the selected region. `log_reg` is the final region log dictionary, which records the parameters, processed protein tar files, kept proteins, and output file paths, etc. used when producing `all_significant_qtls`.


## Contribute

- Issue Tracker: [github.com/nglm/ukbppp-dl/issues](https://github.com/nglm/ukbppp-dl/issues).
- Source Code: [github.com/nglm/ukbppp-dl](github.com/nglm/ukbppp-dl).

## Support

If you are having issues, [please let me know](https://www.uib.no/en/persons/Natacha.Madeleine.Georgette.Galmiche) or [create an issue](https://github.com/nglm/ukbppp-dl/issues).

## License

MIT License. See `LICENSE.txt`.
