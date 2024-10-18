from subprocess import run, PIPE
from tempfile import NamedTemporaryFile

HMMSEARCH_PATH = "hmmsearch"

def hmmsearch(input_file, output_file):
    cmd = f"{HMMSEARCH_PATH} -i {input_file} -f tsv -o {output_file}"
    run(cmd, shell=True, stdout=PIPE, stderr=PIPE)
    return output_file