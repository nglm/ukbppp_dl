import os
from pathlib import Path
from datetime import datetime
import json
from typing import Any

import synapseclient

def save_log(
        log_fname: str,
        log_dict: dict[str, Any],
        overwrite: bool = False,
        add_date: bool = True,
        new_name: bool = False,
        verbose: bool = False,
    ) -> str | None:
    """Save a JSON log file and optionally rename it when a file exists.

    Parameters
    ----------
    log_fname : str
        Target file path for the log file, typically a ``.json`` filename.
    log_dict : dict[str, Any]
        Dictionary of log metadata to serialize. Typical keys include run
        parameters, status flags, and summary values. The ``log_filename`` key
        is added or updated before writing.
    overwrite : bool, default=False
        Whether to overwrite an existing log file at ``log_fname``.
    add_date : bool, default=True
        Whether to append a ``YYYY-MM-DD--HH:MM:SS`` timestamp to the filename
        before saving.
    new_name : bool, default=False
        Whether to create a new timestamped filename when ``log_fname`` already
        exists and ``overwrite`` is ``False``.
    verbose : bool, default=False
        Whether to print status messages during save operations.

    Returns
    -------
    str | None
        The final log filename that was written, or ``None`` if saving was
        skipped.
    """

    if add_date:
        # Create a new filename by adding a suffix to the original filename
        ext = Path(log_fname).suffix
        base = log_fname.split(ext)[0]
        full_date = datetime.today().strftime('%Y-%m-%d--%H:%M:%S')
        log_fname = f"{base}-{full_date}{ext}"

    if os.path.isfile(log_fname):

        # Common message if log exists
        if int(verbose) > 0:
            print(f"Log file {log_fname} already exists.")
        # Added message if log should be overwritten
        if overwrite:
            if int(verbose) > 0:
                print(f"Overwriting log file {log_fname}.")
        # Added message if we then cancel the saving of the log file
        elif not new_name:
            log_fname = None
            if int(verbose) > 0:
                print(f"Not overwriting nor creating new filename. Skipping log saving.")

        # If we create a new name
        else:

            # Create a new filename by adding a suffix to the original filename
            ext = Path(log_fname).suffix
            base = log_fname.split(ext)[0]
            full_date = datetime.today().strftime('%Y-%m-%d--%H:%M:%S')
            log_fname = f"{base}-{full_date}{ext}"

            if int(verbose) > 0:
                print(f"Creating new filename: {log_fname}.")

    if log_fname is not None:

        # Make sure path exists otherwise create it
        p = Path(log_fname)
        p.parent.mkdir(parents=True, exist_ok=True)

        log_dict["log_filename"] = log_fname

        with open(log_fname, 'w') as f_log:
            json.dump(log_dict, f_log, indent=2)

        if int(verbose) > 0:
            print(f"Saved log file to {log_fname}.")

    return log_fname

def delete_files(
        filenames: list[str | None],
        verbose: bool = False,
) -> int:
    """Delete files from disk.

    Parameters
    ----------
    filenames : list[str | None]
        File paths to delete. Typical elements are string paths, with ``None``
        values ignored.
    verbose : bool, default=False
        Whether to print deletion messages.

    Returns
    -------
    int
        Number of files successfully deleted.
    """

    n_deleted = 0
    for fname in filenames:
        if fname is not None and os.path.isfile(fname):
            os.remove(fname)
            n_deleted += 1
            if int(verbose) > 1:
                print(f"Deleting file: {fname}")
    return n_deleted


def download_from_synapse(
        synapse_id: str,
        download_location: str = "./ukb_ppp_dl/data",
        login_kwargs: dict[str, Any] = {},
        verbose: bool = False,
    ) -> tuple[str, bool]:
    """Download a Synapse file unless it already exists locally.

    Parameters
    ----------
    synapse_id : str
        Synapse entity ID for the file to download.
    download_location : str, default="./ukb_ppp_dl/data"
        Directory where the file should be stored. Typical values are local
        project data folders.
    login_kwargs : Dict[str, Any], optional
        Keyword arguments forwarded to
        ``synapseclient.Synapse.login()``. Common keys include
        ``"authToken"`` (``str``), ``"email"`` (``str``) and
        ``"profile"`` (``str``). See `synapseclient.Synapse.login()
        <https://python-docs.synapse.org/en/stable/reference/client/#synapseclient.Synapse.login>`_
        for more information.

        Alternatively, you can configure your Synapse credentials in a
        ``.synapseConfig`` file and leave this argument empty to use the
        default login behaviour. You can find more information about
        ``.synapseConfig`` in the `Synapse documentation
        <https://python-docs.synapse.org/en/stable/tutorials/authentication/>`_.
    verbose : bool, default=False
        Whether to print status messages about login and download decisions.

    Returns
    -------
    tuple[str, bool]
        A tuple containing the expected or downloaded file path and a boolean
        flag indicating whether the download was skipped because the file
        already existed.
    """

    syn = synapseclient.Synapse()
    # Use your synapse account here, as configured in your .synapseConfig file.
    syn.login(**login_kwargs)

    # ----- First check if the file has already been downloaded ------
    # If downloadFile is False, the file will not be downloaded, but we have access to synapse.entity.Entity attributes such as name.
    entity = syn.get(synapse_id, downloadFile=False)

    # Check if the result file already exists in the download location
    expected_file_path = f"{download_location}/{entity.name}"
    if os.path.isfile(expected_file_path):
        skipped = True
        if int(verbose) > 0:
            print(f"File already downloaded: {expected_file_path}. Skipping download.")
    else:
        skipped = False
        # Download the file to the specified location
        entity = syn.get(
            synapse_id, downloadLocation=download_location, downloadFile=True
        )
        expected_file_path = entity.path
        if int(verbose) > 0:
            print(f"Downloaded file: {expected_file_path}")

    return expected_file_path, skipped