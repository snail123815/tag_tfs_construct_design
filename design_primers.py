import argparse
import re
import logging
from pathlib import Path
from tempfile import NamedTemporaryFile

import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq

from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqUtils import MeltingTemp as mt

from modules.find_gene import get_protein
from modules.read_config import read_config
from pyBioinfo_modules.wrappers.hmmer import hmmsearch, read_domtbl
from modules.primer_utils import find_local_best_primer_bind
from pyBioinfo_modules.bio_sequences.bio_features import (
    slice_seq_record_preserve_truncated,
)
from modules.gibson import (
    clean_overhangs_for_gibson,
    make_circular_gibson,
    conjugate_seq_records_gibson,
)

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



def primer_name(prefix, suffix, n):
    return f"{prefix}{n:03d}_{suffix}", n + 1


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
    found_dom_near_term = {}
    target_dom_near_term = {}
    default_term = "N"
    for gene in genes:
        domtbl_gene = domtbl[domtbl["t"] == gene]
        if len(domtbl_gene) == 0:
            logger.warning(f"No hit found for {gene}")
            target_dom_near_term[gene] = default_term
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
        found_dom_near_term[gene] = domain_at_term
        target_dom_near_term[gene] = domain_at_term

        entries.append(domtbl_gene_best_dom[1])

    domtbl = pd.concat(entries, axis=0).reset_index(drop=True)
    domtbl.to_csv(
        args.output / f"{args.name}_domtbl_filtered_best.tsv",
        sep="\t",
        index=False,
    )
    pd.DataFrame(found_dom_near_term, index=["Terminal"]).T.to_csv(
        args.output / f"{args.name}_dom_near_term.tsv", sep="\t", index=True
    )

    # Extract sequences from genome file, prepare for primer3
    cterm_fwd_linker = Seq("GCCAGTATACACTCCGCTAGCG")
    cterm_rev_flag_linker = Seq("TGTCGTGGTCCTTGTAGTCGCCGTCGTGGTCCTTGTAGTC")
    nterm_promoter_fwd_linker = Seq("GCCAGTATACACTCCGCTAGCG")
    nterm_promoter_rev_flag_linker = Seq(
        "TCGTCCTTGTAGTCGATGTCGTGGTCCTTGTAGTCGCCGTCGTGGTCCTTGTAGTC"
    )
    nterm_coding_fwd_flag_linker = Seq("ACCACGACATCGACTACAAGGACGACGACGACAAG")
    nterm_coding_rev_linker = Seq("CAAAGGCCGCTTTTGCGGGATC")

    primers = {}
    primer_pairs = {}
    p_optlength = 20
    wander_range = (-1, 6)
    p_maxlength = p_optlength + wander_range[1]
    promoter_len = 400
    construct_n = configs["construct_id_start"]
    primer_n = configs["oligo_id_start"]
    for gene, hth_term in target_dom_near_term.items():
        construct_id_prefix = f"{configs["construct_prefix"]}{construct_n:03d}"
        construct_n += 1
        gene_rec, contig_i, location = target_sequences[gene]["genes"][0]
        base = contigs[contig_i]
        if hth_term == "N":
            tag_term = "C"
        elif hth_term == "C":
            tag_term = "N"
        else:
            logger.warning(f"Unknown terminal for {gene}")
            continue
        if tag_term == "C":  # HTH in N-terminal, tag C-terminal
            # use the last 20 nucleotides of the gene, reverse complement
            # Check differences around the 20th nucleotide for repeats
            p_rev_anneal, _ = find_local_best_primer_bind(
                gene_rec.seq[-30:-3].reverse_complement()
            )
            if location.strand == 1:
                p_start = location.start - promoter_len
                pext = base[p_start : p_start + p_maxlength]
                ct_product = slice_seq_record_preserve_truncated(
                    base, (p_start, location.end - 3)
                )
            else:
                p_start = location.end + promoter_len
                pext = base[
                    p_start - p_maxlength : p_start
                ].reverse_complement()
                ct_product = slice_seq_record_preserve_truncated(
                    base, (location.start + 3, p_start)
                ).reverse_complement(
                    id=f"{base.id}_rev_comp",
                    name=f"{base.name}_rev_comp",
                    description=base.description,
                    annotations={
                        "molecule_type": base.annotations["molecule_type"]
                    },
                )
            p_fwd_anneal, _ = find_local_best_primer_bind(
                pext, p_optlength, wander_range
            )
            pf_name, primer_n = primer_name(
                configs["oligo_prefix"], f"{gene}_Ct_fwd", primer_n
            )
            pr_name, primer_n = primer_name(
                configs["oligo_prefix"], f"{gene}_Ct_rev", primer_n
            )
            pf = cterm_fwd_linker + p_fwd_anneal
            pr = cterm_rev_flag_linker + p_rev_anneal
            ct_product = (
                cterm_fwd_linker
                + ct_product
                + cterm_rev_flag_linker.reverse_complement()
            )
            primers[pf_name] = {
                "anneal": p_fwd_anneal,
                "linker": cterm_fwd_linker,
                "with_linker": pf,
                "full_len": len(pf),
                "anneal_len": len(p_fwd_anneal),
            }
            primers[pr_name] = {
                "anneal": p_rev_anneal,
                "linker": cterm_rev_flag_linker,
                "with_linker": pr,
                "full_len": len(pr),
                "anneal_len": len(p_rev_anneal),
            }
            primer_pairs[f"{construct_id_prefix}_{gene}_Ct"] = {
                "fwd": pf_name,
                "fwd_TM": f"{mt.Tm_NN(p_fwd_anneal):.2f}",
                "rev": pr_name,
                "rev_TM": f"{mt.Tm_NN(p_rev_anneal):.2f}",
                "product_size": (
                    (location.end - location.start)
                    - 3
                    + promoter_len
                    + len(cterm_fwd_linker)
                    + len(cterm_rev_flag_linker)
                ),
                "fwd_seq": pf,
                "fwd_full_len": len(pf),
                "rev_seq": pr,
                "rev_full_len": len(pr),
                "fwd_anneal": p_fwd_anneal,
                "fwd_anneal_len": len(p_fwd_anneal),
                "rev_anneal": p_rev_anneal,
                "rev_anneal_len": len(p_rev_anneal),
                "product": ct_product,
            }
        elif tag_term == "N":  # HTH in C-terminal, tag N-terminal
            if location.strand == 1:
                promoter_start = location.start + 3 - promoter_len
                promoter_end = location.start + 3
                coding_start = location.start + 3
                coding_end = location.end
                p_promoter_fwd_anneal, _ = find_local_best_primer_bind(
                    base[promoter_start : promoter_start + p_maxlength],
                )
                p_promoter_rev_anneal, _ = find_local_best_primer_bind(
                    base[
                        promoter_end - p_maxlength : promoter_end
                    ].reverse_complement(),
                )
                p_coding_fwd_anneal, _ = find_local_best_primer_bind(
                    base[coding_start : coding_start + p_maxlength],
                )
                p_coding_rev_anneal, _ = find_local_best_primer_bind(
                    base[
                        coding_end - p_maxlength : coding_end
                    ].reverse_complement(),
                )
                nt_promoter_product_no_linker = (
                    slice_seq_record_preserve_truncated(
                        base, (promoter_start, promoter_end)
                    )
                )
                nt_coding_product_no_linker = (
                    slice_seq_record_preserve_truncated(
                        base, (coding_start, coding_end)
                    )
                )
            else:
                promoter_start = location.end - 3 + promoter_len
                promoter_end = location.end - 3
                coding_start = location.end - 3
                coding_end = location.start
                p_promoter_fwd_anneal, _ = find_local_best_primer_bind(
                    base[
                        promoter_start - p_maxlength : promoter_start
                    ].reverse_complement(),
                )
                p_promoter_rev_anneal, _ = find_local_best_primer_bind(
                    base[promoter_end : promoter_end + p_maxlength],
                )
                p_coding_fwd_anneal, _ = find_local_best_primer_bind(
                    base[
                        coding_start - p_maxlength : coding_start
                    ].reverse_complement(),
                )
                p_coding_rev_anneal, _ = find_local_best_primer_bind(
                    base[coding_end : coding_end + p_maxlength],
                )
                nt_promoter_product_no_linker = (
                    slice_seq_record_preserve_truncated(
                        base, (promoter_end, promoter_start)
                    ).reverse_complement(
                        id=f"{base.id}_rev_comp",
                        name=f"{base.name}_rev_comp",
                        description=base.description,
                        annotations={
                            "molecule_type": base.annotations["molecule_type"]
                        },
                    )
                )
                nt_coding_product_no_linker = (
                    slice_seq_record_preserve_truncated(
                        base, (coding_end, coding_start)
                    ).reverse_complement(
                        id=f"{base.id}_rev_comp",
                        name=f"{base.name}_rev_comp",
                        description=base.description,
                        annotations={
                            "molecule_type": base.annotations["molecule_type"]
                        },
                    )
                )

            nt_promoter_product = (
                nterm_promoter_fwd_linker
                + nt_promoter_product_no_linker
                + nterm_promoter_rev_flag_linker.reverse_complement()
            )
            nt_coding_product = (
                nterm_coding_fwd_flag_linker
                + nt_coding_product_no_linker
                + nterm_coding_rev_linker.reverse_complement()
            )
            ppf_name, primer_n = primer_name(
                configs["oligo_prefix"], f"{gene}_Nt_promoter_fwd", primer_n
            )
            ppr_name, primer_n = primer_name(
                configs["oligo_prefix"], f"{gene}_Nt_promoter_rev", primer_n
            )
            ppf = nterm_promoter_fwd_linker + p_promoter_fwd_anneal
            ppr = nterm_promoter_rev_flag_linker + p_promoter_rev_anneal
            primers[ppf_name] = {
                "anneal": p_promoter_fwd_anneal,
                "linker": nterm_promoter_fwd_linker,
                "with_linker": ppf,
                "full_len": len(ppf),
                "anneal_len": len(p_promoter_fwd_anneal),
            }
            primers[ppr_name] = {
                "anneal": p_promoter_rev_anneal,
                "linker": nterm_promoter_rev_flag_linker,
                "with_linker": ppr,
                "full_len": len(ppr),
                "anneal_len": len(p_promoter_rev_anneal),
            }
            primer_pairs[f"{construct_id_prefix}_{gene}_Nt_promoter"] = {
                "fwd": ppf_name,
                "fwd_TM": f"{mt.Tm_NN(p_promoter_fwd_anneal):.2f}",
                "rev": ppr_name,
                "rev_TM": f"{mt.Tm_NN(p_promoter_rev_anneal):.2f}",
                "product_size": (
                    abs(promoter_start - promoter_end)
                    + len(nterm_promoter_fwd_linker)
                    + len(nterm_promoter_rev_flag_linker)
                ),
                "fwd_seq": ppf,
                "fwd_full_len": len(ppf),
                "rev_seq": ppr,
                "rev_full_len": len(ppr),
                "fwd_anneal": p_promoter_fwd_anneal,
                "fwd_anneal_len": len(p_promoter_fwd_anneal),
                "rev_anneal": p_promoter_rev_anneal,
                "rev_anneal_len": len(p_promoter_rev_anneal),
                "product": nt_promoter_product,
            }

            pcf_name, primer_n = primer_name(
                configs["oligo_prefix"], f"{gene}_Nt_coding_fwd", primer_n
            )
            pcr_name, primer_n = primer_name(
                configs["oligo_prefix"], f"{gene}_Nt_coding_rev", primer_n
            )
            pcf = nterm_coding_fwd_flag_linker + p_coding_fwd_anneal
            pcr = nterm_coding_rev_linker + p_coding_rev_anneal
            primers[pcf_name] = {
                "anneal": p_coding_fwd_anneal,
                "linker": nterm_coding_fwd_flag_linker,
                "with_linker": pcf,
                "full_len": len(pcf),
                "anneal_len": len(p_coding_fwd_anneal),
            }
            primers[pcr_name] = {
                "anneal": p_coding_rev_anneal,
                "linker": nterm_coding_rev_linker,
                "with_linker": pcr,
                "full_len": len(pcr),
                "anneal_len": len(p_coding_rev_anneal),
            }
            primer_pairs[f"{construct_id_prefix}_{gene}_Nt_coding"] = {
                "fwd": pcf_name,
                "fwd_TM": f"{mt.Tm_NN(p_coding_fwd_anneal):.2f}",
                "rev": pcr_name,
                "rev_TM": f"{mt.Tm_NN(p_coding_rev_anneal):.2f}",
                "product_size": (
                    abs(coding_start - coding_end)
                    + len(nterm_coding_fwd_flag_linker)
                    + len(nterm_coding_rev_linker)
                ),
                "fwd_seq": pcf,
                "fwd_full_len": len(pcf),
                "rev_seq": pcr,
                "rev_full_len": len(pcr),
                "fwd_anneal": p_coding_fwd_anneal,
                "fwd_anneal_len": len(p_coding_fwd_anneal),
                "rev_anneal": p_coding_rev_anneal,
                "rev_anneal_len": len(p_coding_rev_anneal),
                "product": nt_coding_product,
            }
        else:
            logger.warning(f"Unknown terminal for {gene}")

    # Write primers to file
    pd.DataFrame(primers).T.to_csv(
        args.output / f"{args.name}_primers.tsv", sep="\t", index=True
    )
    primer_pairs_df = pd.DataFrame(primer_pairs).T
    primer_pairs_df["product"] = primer_pairs_df["product"].apply(
        lambda x: str(x.seq)
    )
    primer_pairs_df.index.name = "fragment"
    primer_pairs_df.to_csv(
        args.output / f"{args.name}_pairs.tsv", sep="\t", index=True
    )

    # Make constructs
    # 1. remove BsaI overhangs
    ct_backbone = clean_overhangs_for_gibson(
        SeqIO.read(configs["backbone_cterm"], "genbank")
    )
    nt_backbone = clean_overhangs_for_gibson(
        SeqIO.read(configs["backbone_nterm"], "genbank")
    )
    # 2. find overlaps and remove from backbone them
    # 3. concatenate sequence, output as genbank
    for gene, hth_term in target_dom_near_term.items():
        if hth_term == "N":
            tag_term = "C"
        elif hth_term == "C":
            tag_term = "N"
        else:
            logger.warning(f"Unknown terminal for {gene}")
            continue
        primers_per_gene = {}
        for p in primers:
            if gene in p:
                primers_per_gene[p] = primers[p]
        frags = {}
        for pp in primer_pairs:
            if gene in pp:
                frags[pp] = primer_pairs[pp]["product"]
        if tag_term == "C":
            assert (
                len(frags) == 1
            ), "Multiple fragments found for C trem tagging"
            construct_name = [p for p in frags][0]
            frag = frags[construct_name]
            construct = make_circular_gibson(
                conjugate_seq_records_gibson(frag, ct_backbone)
            )
            construct.id = construct_name
        elif tag_term == "N":
            assert len(frags) == 2, "Needs 2 fragments N trem tagging"
            promoter_name = [p for p in frags if "_promoter" in p][0]
            coding_name = [p for p in frags if "_coding" in p][0]
            frag_promoter = frags[promoter_name]
            frag_coding = frags[coding_name]
            construct = make_circular_gibson(
                conjugate_seq_records_gibson(
                    conjugate_seq_records_gibson(frag_promoter, frag_coding),
                    nt_backbone,
                )
            )
            construct.id = re.sub(r"_promoter$", "", promoter_name)
        else:
            logger.warning(f"Unknown terminal for {gene}")
        construct.name = construct.id
        # Add primers to construct
        for f in construct.features:
            if any(
                [f.qualifiers[q][0].startswith("ori") for q in f.qualifiers]
            ):
                if f.location.strand == -1:
                    zero_location = f.location.end
                    strand = -1
                else:
                    zero_location = f.location.start
                    strand = 1
        a = slice_seq_record_preserve_truncated(construct, (0, zero_location))
        b = slice_seq_record_preserve_truncated(
            construct, (zero_location, len(construct))
        )
        construct_meta = {
            "id": construct.id,
            "name": construct.name,
            "description": construct.description,
            "annotations": construct.annotations,
        }
        construct = b + a
        if strand == -1:
            construct = construct.reverse_complement()
        construct.id = construct_meta["id"]
        construct.name = construct_meta["name"]
        construct.description = construct_meta["description"]
        construct.annotations = construct_meta["annotations"]
        for p in primers_per_gene:
            # find location of primer in construct
            pseq = primers_per_gene[p]["with_linker"]
            strand = 1
            ploc = construct.seq.find(pseq)
            if ploc == -1:
                ploc = construct.seq.find(pseq.reverse_complement())
                strand = -1
            if ploc == -1:
                logger.warning(f"Primer {p} not found in construct")
                continue
            primer_feat = SeqFeature(
                FeatureLocation(ploc, ploc + len(pseq), strand=strand),
                type="primer_bind",
                qualifiers={
                    "label": [p],
                    "primer": [pseq],
                    "overhang": [primers[p]["linker"]],
                    "anneal": [primers[p]["anneal"]],
                },
            )
            construct.features.append(primer_feat)
        construct.features = sorted(
            construct.features, key=lambda x: x.location.start
        )
        SeqIO.write(construct, f"{args.output/construct.id}.gbk", "genbank")


if __name__ == "__main__":
    main()
