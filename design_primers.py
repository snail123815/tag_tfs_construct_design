import argparse
import json
import logging
import numpy as np
from pathlib import Path
from tempfile import NamedTemporaryFile

from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd

from modules.find_gene import get_protein
from modules.read_config import read_config
from pyBioinfo_modules.wrappers.hmmer import hmmsearch, read_domtbl

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
    hmmsearch(configs["hmmsearch"], protein_faa.name, args.hmm, hmmtblout)

    # Parse hmmtblout
    domtbl = read_domtbl(hmmtblout, t_e=configs["gather_threshold_e"])

    # Sort by E value and find the best hit
    entries = []
    target_dom_near_term = {}
    for gene in genes:
        domtbl_gene = domtbl[domtbl["t"] == gene]
        if len(domtbl_gene) == 0:
            logger.warning(f"No hit found for {gene}")
            continue
        domtbl_gene_doms = []
        # Make "HTH" domain the priority
        if any(["HTH" in d for d in domtbl_gene["q"]]):
            found_doms = domtbl_gene[domtbl_gene["q"].str.contains("HTH")][
                "q"
            ].unique()
        else:
            found_doms = domtbl_gene["q"].unique()
        # Find the best domain
        for dom in found_doms:
            domtbl_dom = domtbl_gene[domtbl_gene["q"] == dom]
            dom_min_e = domtbl_dom["dom_i_E"].min()
            dom_e_transformed = -np.log10(domtbl_dom["dom_i_E"])
            # Only include domains with E value higher than 10 magnitudes of
            # the best hit
            domtbl_dom = domtbl_dom[
                dom_e_transformed >= (dom_e_transformed.max() / 10)
            ]
            dom_start = domtbl_dom["ali_from"].min()
            dom_end = domtbl_dom["ali_to"].max()
            domtbl_gene_doms.append(
                [dom, domtbl_dom, dom_min_e, dom_start, dom_end]
            )
        domtbl_gene_doms = sorted(domtbl_gene_doms, key=lambda x: x[2])
        domtbl_gene_best_dom = domtbl_gene_doms[0]
        start = domtbl_gene_best_dom[3]
        end = domtbl_gene_best_dom[4]
        ndist = start - 1
        cdist = domtbl_gene_best_dom[1]["tlen"].iloc[0] - end
        logging.info(f"Gene: {gene}")
        logging.info(
            f"C-terminal distance: {ndist}, N-terminal distance: {cdist}"
        )
        domain_at_term = "N" if ndist <= cdist else "C"
        logging.info(f"Domain at terminal: {domain_at_term}")
        target_dom_near_term[gene] = domain_at_term

        entries.append(domtbl_gene_best_dom[1])

    domtbl = pd.concat(entries, axis=0).reset_index(drop=True)
    domtbl.to_csv(
        args.output / f"{args.name}_domtbl_filtered_best.tsv",
        sep="\t",
        index=False,
    )
    pd.DataFrame(target_dom_near_term, index=["Terminal"]).T.to_csv(
        args.output / f"{args.name}_dom_near_term.tsv", sep="\t", index=True
    )

    # Extract sequences from genome file, prepare for primer3
    nterm_flag_linker = Seq("GTCCTTGTAGTCGCCGTCGTGGTCCTTGTAGTC")
    primers_tag_cterm = {}
    primers_tag_nterm = {}
    for gene, term in target_dom_near_term.items():
        gene_rec, contig_i, location = target_sequences[gene]["genes"][0]  #
        base = contigs[contig_i]
        if term == "N":
            # use the last 20 nucleotides of the gene, reverse complement
            p_rev_anneal = gene_rec.seq[-30:-3].reverse_complement()
            # Check differences around the 20th nucleotide for repeats
            found_diff = False
            for i in list(range(7)):
                if p_rev_anneal[-5 + i] != p_rev_anneal[-5 + i + 1]:
                    found_diff = True
                    p_rev_anneal = p_rev_anneal[: -5 + i + 1]
                    break
            if not found_diff:
                p_rev_anneal = p_rev_anneal[:-5]

            target_location_floc = (
                location.start - 400
                if location.strand == 1
                else location.end + 400
            )
            p_fwd_anneal = (
                base[target_location_floc : target_location_floc + 20]
                if location.strand == 1
                else base[
                    target_location_floc - 20 : target_location_floc
                ].reverse_complement()
            )
            primers_tag_cterm[gene] = {
                "fwd_anneal": p_fwd_anneal,
                "rev_anneal": p_rev_anneal,
            }
        else:
            p_promoter_fwd_anneal = (
                base[location.start - 400 : location.start - 380]
                if location.strand == 1
                else base[
                    location.end + 380 : location.end + 400
                ].reverse_complement()
            )
            p_promoter_rev_anneal = (
                base[location.start - 20 : location.start].reverse_complement()
                if location.strand == 1
                else base[location.end : location.end + 20]
            )
            p_gene_fwd_anneal = (
                base[location.start + 3 : location.start + 23]
                if location.strand == 1
                else base[
                    location.end - 23 : location.end - 3
                ].reverse_complement()
            )
            p_gene_rev_anneal = (
                base[location.end - 20 : location.end].reverse_complement()
                if location.strand == 1
                else base[location.start : location.start + 20]
            )

            primers_tag_nterm[gene] = {
                "promoter_fwd_anneal": p_promoter_fwd_anneal,
                "promoter_rev_anneal": p_promoter_rev_anneal,
                "gene_fwd_anneal": p_gene_fwd_anneal,
                "gene_rev_anneal": p_gene_rev_anneal,
            }
            
    # primer3(data["genes"], args.output, configs)


if __name__ == "__main__":
    main()
