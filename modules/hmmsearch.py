import subprocess
from tempfile import NamedTemporaryFile
from pathlib import Path
import logging

logger = logging.getLogger(__name__)


def hmmsearch(
    HMMSEARCH_CMD: Path, protein_faa: Path, hmm_file: Path, output_file: Path
):
    cmd = [
        str(HMMSEARCH_CMD),
        "--domtblout",
        str(output_file),
        str(hmm_file),
        str(protein_faa),
    ]
    search_run = subprocess.run(cmd, capture_output=True)
    try:
        search_run.check_returncode()
    except subprocess.CalledProcessError as e:
        logger.error(
            f"Error running hmmsearch: {search_run.stderr.decode()}"
        )
        raise e
    return output_file
