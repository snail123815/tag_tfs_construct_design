from copy import deepcopy
from math import ceil

from Bio.SeqFeature import (
    AfterPosition,
    BeforePosition,
    FeatureLocation,
    SeqFeature,
)
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


def add_seq_to_SeqRecord_as_feature(
    seq_record: SeqRecord,
    target_seq: Seq,
    feature_type: str = "feature",
    qualifiers: dict = {},
) -> SeqRecord:
    """
    Add sequence to the end of the sequence record. Automatically finds the
    correct strand and location for the target sequence.
    If a primer bind sequence is provided, it will check if there is
    5' overhang, and avoid that.
    """

    if (
        feature_type == "primer_bind" and "overhang" in qualifiers
    ):  # 5' overhang
        try:
            overhang = int(qualifiers["overhang"][0])
            overhang = Seq("N" * overhang)
        except ValueError:
            overhang = Seq(qualifiers["overhang"][0])
            assert (
                overhang in target_seq
            ), "Overhang not found in target sequence"
        anneal_seq = target_seq[len(overhang) :]
    else:
        anneal_seq = target_seq

    seq_seq = seq_record.seq
    seq_len = len(seq_record)
    anneal_len = len(anneal_seq)
    on_strand = 1
    if anneal_seq.reverse_complement() in seq_seq:
        # located on -1 strand
        on_strand = -1
    else:
        raise ValueError("Oligo not found in sequence")
    match_seq = (
        anneal_seq if on_strand == 1 else anneal_seq.reverse_complement()
    )
    for i in range(seq_len - anneal_len + 1):
        if seq_seq[i : i + anneal_len] == match_seq:
            location = FeatureLocation(i, i + anneal_len, strand=on_strand)
            feature = SeqFeature(
                location=location,
                type=feature_type,
                qualifiers=qualifiers,
            )
            seq_record.features.append(feature)
            return seq_record
    raise ValueError("Oligo not found in sequence")


def slice_seq_record_preserve_truncated(
    seq_record: SeqRecord, slice_tuple: tuple[int, int]
) -> SeqRecord:
    sliced_rec = seq_record[slice_tuple[0] : slice_tuple[1]]
    truncated_features_left = []
    truncated_features_right = []
    for feature in seq_record.features:
        if feature.location.start < slice_tuple[0] < feature.location.end:
            feat = deepcopy(feature)
            feat.location = (
                FeatureLocation(
                    BeforePosition(slice_tuple[0]),
                    (
                        feat.location.end
                        if feat.location.end <= slice_tuple[1]
                        else AfterPosition(slice_tuple[1])
                    ),
                    feat.location.strand,
                )
                - slice_tuple[0]
            )
            feat.id = f"{feature.id}_truncated"
            feat.qualifiers["truncated"] = ["left"]
            truncated_features_left.append(feat)
            continue
        if feature.location.start < slice_tuple[1] < feature.location.end:
            feat = deepcopy(feature)
            feat.location = (
                FeatureLocation(
                    (
                        feat.location.start
                        if feat.location.start >= slice_tuple[0]
                        else BeforePosition(slice_tuple[0])
                    ),
                    AfterPosition(slice_tuple[1]),
                    feat.location.strand,
                )
                - slice_tuple[0]
            )
            feat.id = f"{feature.id}_truncated"
            feat.qualifiers["truncated"] = ["right"]
            truncated_features_right.append(feat)
    sliced_rec.features = (
        truncated_features_left + sliced_rec.features + truncated_features_right
    )
    return sliced_rec


def getSpanFetures(genome, startOri, endOri, expand=20000):
    newStart = max(0, startOri - expand)
    newEnd = min(len(genome), endOri + expand)
    sourceSeq = genome[newStart:newEnd]

    start = startOri - newStart
    end = start + (endOri - startOri)

    spanFeats = []
    for feat in sourceSeq.features:
        spanStart = start in feat
        spanEnd = end in feat
        # feat.__contains__(self, value)
        # Check if an integer position is within the feature.

        spanFeat = SeqFeature(type=feat.type)

        if spanStart and spanEnd:
            # Target position is inside feature
            if feat.type == "CDS":
                # calculate correct start and end location to make it inframe
                newStart = (3 - abs(start - feat.location.start)) % 3
                newEnd = end - start - abs(end - feat.location.start) % 3
            else:
                newStart = 0
                newEnd = end - start
            spanFeat.location = FeatureLocation(
                newStart, newEnd, strand=feat.location.strand
            )
            for key in feat.qualifiers:
                if key in [
                    "gene_synonym",
                    "locus_tag",
                    "product",
                    "protein_id",
                    "db_xref",
                    "mol_type",
                ]:
                    spanFeat.qualifiers[key] = [
                        f"{keyStr} (slice)" for keyStr in feat.qualifiers[key]
                    ]
                elif key == "translation":
                    spanFeat.qualifiers[key] = list(feat.qualifiers[key])
                    if feat.location.strand == 1:
                        cutPointA = ceil((start - feat.location.start) / 3)
                        cutPointB = (end - feat.location.start) // 3
                    else:
                        len(spanFeat.qualifiers[key][0])
                        cutPointA = (
                            len(spanFeat.qualifiers[key][0])
                            + 1
                            - (end - feat.location.start) // 3
                        )
                        cutPointB = (
                            len(spanFeat.qualifiers[key][0])
                            + 1
                            - ceil((start - feat.location.start) / 3)
                        )
                    spanFeat.qualifiers[key][0] = spanFeat.qualifiers[key][0][
                        cutPointA:cutPointB
                    ]
                else:
                    spanFeat.qualifiers[key] = feat.qualifiers[key]
            spanFeats.append(spanFeat)

        elif spanStart:
            # Start position inside feature, feature ends in this range
            if feat.type == "CDS":
                newStart = (3 - abs(start - feat.location.start)) % 3
            else:
                newStart = 0
            newEnd = feat.location.end - start
            spanFeat.location = FeatureLocation(
                newStart, newEnd, strand=feat.location.strand
            )
            for key in feat.qualifiers:
                if key in [
                    "gene_synonym",
                    "locus_tag",
                    "product",
                    "protein_id",
                    "db_xref",
                    "mol_type",
                ]:
                    spanFeat.qualifiers[key] = [
                        f"{keyStr} (right part)"
                        for keyStr in feat.qualifiers[key]
                    ]
                elif key == "translation":
                    spanFeat.qualifiers[key] = [
                        keyStr for keyStr in feat.qualifiers[key]
                    ]
                    if feat.location.strand == 1:
                        cutPoint = (
                            len(spanFeat.qualifiers[key][0])
                            - ceil(len(spanFeat) / 3)
                            + 1
                        )
                        spanFeat.qualifiers[key][0] = spanFeat.qualifiers[key][
                            0
                        ][cutPoint:]
                    else:
                        cutPoint = len(spanFeat) // 3
                        spanFeat.qualifiers[key][0] = spanFeat.qualifiers[key][
                            0
                        ][:cutPoint]
                else:
                    spanFeat.qualifiers[key] = feat.qualifiers[key]
            spanFeats.append(spanFeat)

        elif spanEnd:
            # End position inside feature, feature ends in this range
            if feat.type == "CDS":
                newEnd = end - start - abs(end - feat.location.start) % 3
            else:
                newEnd = end - start
            newStart = feat.location.start - start
            spanFeat.location = FeatureLocation(
                newStart, newEnd, strand=feat.location.strand
            )
            for key in feat.qualifiers:
                if key in [
                    "gene_synonym",
                    "locus_tag",
                    "product",
                    "protein_id",
                    "db_xref",
                    "mol_type",
                ]:
                    spanFeat.qualifiers[key] = [
                        f"{keyStr} (left part)"
                        for keyStr in feat.qualifiers[key]
                    ]
                elif key == "translation":
                    spanFeat.qualifiers[key] = list(feat.qualifiers[key])
                    if feat.location.strand == 1:
                        cutPoint = len(spanFeat) // 3
                        spanFeat.qualifiers[key][0] = spanFeat.qualifiers[key][
                            0
                        ][:cutPoint]
                    else:
                        cutPoint = ceil((len(feat) - len(spanFeat)) / 3)
                        spanFeat.qualifiers[key][0] = spanFeat.qualifiers[key][
                            0
                        ][cutPoint:]
                else:
                    spanFeat.qualifiers[key] = feat.qualifiers[key]
            spanFeats.append(spanFeat)

        else:
            # Not in range, ignore
            continue

    return spanFeats


# getSpanFetures
