from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

import logging

logger = logging.getLogger(__name__)


ID_QUALIFIERS = [
    "locus_tag",
    "gene",
    "protein_id",
    "gene_synonym",
    "product",
]


def get_dna(contig, location):
    """
    Get the DNA sequence of a location in a contig.
    Will return the reverse complement if strand = -1.
    """
    seq = contig.seq[location.start : location.end]
    if location.strand == -1:
        seq = seq.reverse_complement()
    return seq


def get_protein(genome, genes, codon_table=11):
    target_sequences = {}
    # Read the genome file and find targets
    contigs = list(SeqIO.parse(genome, "genbank"))
    for contig_i, contig in enumerate(contigs):
        for feature in contig.features:
            if feature.type in ["CDS", "gene"]:
                id = None
                for qualifier in ID_QUALIFIERS:
                    qid = feature.qualifiers.get(qualifier, [""])[0]
                    if qid in genes:
                        id = qid
                        break
                if id is None:
                    continue
                if feature.type == "gene":
                    gene = SeqRecord(
                        get_dna(contig, feature.location),
                        id=id,
                        description="",
                    )
                    target_sequences.setdefault(id, {}).setdefault(
                        "genes", []
                    ).append((gene, contig_i, feature.location))
                else:
                    faa = feature.qualifiers.get("translation", [""])[0]
                    if faa == "":
                        fna = contig.seq[
                            feature.location.start : feature.location.end
                        ]
                        if feature.strand == -1:
                            fna = fna.reverse_complement()
                        faa = fna.translate(table=codon_table, to_stop=True)
                        if faa * 3 != len(fna):
                            logger.warning(
                                f"Translation of gene {id} does not match "
                                "gene length."
                            )
                    faa = Seq(faa)
                    protein = SeqRecord(
                        faa,
                        id=id,
                        description="",
                    )
                    target_sequences.setdefault(id, {}).setdefault(
                        "proteins", []
                    ).append((protein, contig_i, feature.location))

    # Check if all genes were found, if some are found multiple times, etc.
    for gene in genes:
        no_genes = False
        no_proteins = False
        if gene not in target_sequences:
            logger.warning(f"Gene {gene} not found in genome.")
            continue
        if "genes" not in target_sequences[gene]:
            logger.warning(f"Gene {gene} not found in genome.")
            no_genes = True
        if "proteins" not in target_sequences[gene]:
            logger.warning(f"CDS of {gene} not found in genome.")
            no_proteins = True
        if len(target_sequences[gene]["genes"]) > 1:
            logger.warning(f"Gene {gene} found multiple times in genome.")
        if len(target_sequences[gene]["proteins"]) > 1:
            logger.warning(f"CDS of {gene} found multiple times in genome.")
        if no_genes:
            for protein, contig_i, location in target_sequences[gene][
                "proteins"
            ]:
                target_sequences[gene].setdefault("genes", []).append(
                    (
                        SeqRecord(
                            get_dna(contigs[contig_i], location),
                            id=id,
                            description="",
                        ),
                        contig_i,
                        location,
                    )
                )
        if no_proteins:
            for gene, contig_i, location in target_sequences[gene]["genes"]:
                target_sequences[gene].setdefault("proteins", []).append(
                    (
                        SeqRecord(
                            get_dna(contigs[contig_i], location).translate(
                                table=codon_table, to_stop=True
                            ),
                            id=id,
                            description="",
                        ),
                        contig_i,
                        location,
                    )
                )

    return contigs, target_sequences
