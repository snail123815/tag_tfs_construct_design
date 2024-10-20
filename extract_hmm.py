import argparse
from pathlib import Path
import logging

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
)


def argparser():
    parser = argparse.ArgumentParser(
        description="Extract listed HMM profiles from Pfam-A.hmm"
    )
    parser.add_argument(
        "--pfam", type=Path, required=True, help="Path to Pfam-A.hmm"
    )
    parser.add_argument(
        "--accessions",
        type=str,
        required=True,
        help="Path of accession file, or a string of accessions separated by ','",
    )
    parser.add_argument(
        "--output", type=Path, required=True, help="Path to output file"
    )

    return parser

def extract_hmm_profiles(input_file: Path, accessions: set, output_file: Path):
    extracted_accessions = {}
    with input_file.open("r") as infile, output_file.open("w") as outfile:
        block = []
        acc_ver = ""
        acc = ""
        is_target = False
        for line in infile:
            block.append(line)
            if line.startswith("ACC"):
                acc_ver = line.split()[1]
                acc = acc_ver.split(".")[0]
                if acc in accessions:
                    is_target = True
            if line.startswith("//"):
                if is_target:
                    for l in block:
                        outfile.write(l)
                    extracted_accessions.setdefault(acc, []).append(acc_ver)
                acc_ver, acc, is_target = "", "", False
                block = []
    if len(accessions) > len(extracted_accessions):
        not_found = set(accessions) - set(extracted_accessions)
        logging.warning(
            "Some accessions were not found in the HMM file.\n"
            f"Accessions not found: {', '.join(not_found)}"
        )
    for acc in extracted_accessions:
        if len(extracted_accessions[acc]) > 1:
            logging.info(
                f"Accession {acc} found {len(extracted_accessions[acc])} times, "
                "all are written to the output file."
            )
    return extracted_accessions

def main():
    parser = argparser()
    args = parser.parse_args()
    assert args.pfam.is_file(), "Pfam-A.hmm file not found"
    assert args.output.parent.is_dir(), "Output folder not found"
    if Path(args.accessions).is_file():
        accessions = set()
        with open(args.accessions, "r") as infile:
            for line in infile:
                if line.startswith("Accession"):
                    continue
                accessions.add(line.split("\t")[0].strip())
    else:
        accessions = set(args.accessions.split(","))

    extract_hmm_profiles(
        args.pfam, accessions, args.output
    )


if __name__ == "__main__":
    main()
