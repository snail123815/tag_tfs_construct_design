import argparse
import json
import logging
from pathlib import Path
from tempfile import NamedTemporaryFile

from Bio import SeqIO

from modules.find_gene import get_protein
from modules.read_config import read_config
from pyBioinfo_modules.wrappers.hmmer import hmmsearch

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
)
logger = logging.getLogger(__name__)


def argparser():
    parser = argparse.ArgumentParser(
        description="Design primers for given genes in a genome."
    )
    parser.add_argument(
        "--genome", type=Path, required=True, help="Path of genome genbank file"
    )
    parser.add_argument(
        "--genes",
        type=Path,
        required=True,
        help='Path to a file of list of genes, or a string of genes, separated by ","',
    )
    parser.add_argument(
        "--hmm", type=Path, required=True, help="Path to hmm file"
    )
    parser.add_argument(
        "--name", type=str, default="primer_design", help="Name of the output"
    )
    parser.add_argument(
        "--output", type=Path, required=True, help="Path to output folder"
    )
    parser.add_argument(
        "--config", type=Path, required=True, help="Path to config file"
    )

    return parser


def main():
    parser = argparser()
    args = parser.parse_args()
    assert args.genome.is_file(), "Genome file not found"
    assert args.config.is_file(), "Config file not found"
    if args.genes.is_file():
        with open(args.genes) as f:
            genes = f.read().splitlines()
    else:
        genes = str(args.genes).split(",")

    configs = read_config(args.config)

    contigs, target_sequences = get_protein(args.genome, genes)

    for gene, data in target_sequences.items():
        if len(data["genes"]) > 1:
            logger.warning(
                f"Multiple genes found for {gene}, using the first one"
            )
        if len(data["proteins"]) > 1:
            logger.warning(
                f"Multiple proteins found for {gene}, using the first one"
            )

    # Write protein sequences as temporary fasta file (.faa)
    proteins = []
    for gene, data in target_sequences.items():
        proteins.append(data["proteins"][0][0])
    protein_faa = NamedTemporaryFile(mode="w", delete=False)
    SeqIO.write(proteins, protein_faa.name, "fasta")

    # Run hmmsearch
    args.output.mkdir(exist_ok=True)
    hmmtblout = args.output / f"{args.name}_{args.hmm.stem}.hmmtblout"
    hmmsearch_run = hmmsearch(
        configs["hmmsearch"], protein_faa.name, args.hmm, hmmtblout
    )

    # primer3(data["genes"], args.output, configs)


if __name__ == "__main__":
    main()
