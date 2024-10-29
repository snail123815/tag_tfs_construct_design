from datetime import datetime
import logging
from Bio.SeqFeature import SeqFeature, FeatureLocation, ExactPosition
from Bio.SeqRecord import SeqRecord
from pyBioinfo_modules.bio_sequences.bio_features import slice_seq_record_preserve_truncated

logger = logging.getLogger(__name__)

def clean_overhangs_for_gibson(seq_record: SeqRecord) -> SeqRecord:
    first_feature = seq_record.features[0]
    last_feature = seq_record.features[-1]
    if "overhang" in first_feature.qualifiers["note"][0]:
        if first_feature.location.strand == 1:
            seq_record = slice_seq_record_preserve_truncated(
                seq_record, (len(first_feature), len(seq_record))
            )
    if "overhang" in last_feature.qualifiers["note"][0]:
        if last_feature.location.strand == -1:
            seq_record = slice_seq_record_preserve_truncated(
                seq_record, (0, len(seq_record) - len(last_feature))
            )
    return seq_record

def make_circular_gibson(seq_record: SeqRecord) -> SeqRecord:
    for i in range(15, int(len(seq_record) / 2)):
        if seq_record.seq[-i:] == seq_record.seq[:i]:
            break
    seq_record = slice_seq_record_preserve_truncated(
        seq_record, (i, len(seq_record))
    )
    linker_feature = SeqFeature(
        FeatureLocation(
            ExactPosition(len(seq_record) - i),
            ExactPosition(len(seq_record)),
            strand=None,
        ),
        type="misc_feature",
        qualifiers={
            "note": ["Gibson assembly linker"],
            "standard_name": ["Ligation"],
        },
    )
    seq_record.features.append(linker_feature)
    seq_record.annotations["topology"] = "circular"
    return seq_record

def conjugate_seq_records_gibson(
    left_record: SeqRecord, right_record: SeqRecord
) -> SeqRecord:
    left_seq = left_record.seq
    right_seq = right_record.seq
    annotations = {
        "topology": "linear",
        "molecule_type": "DNA",
        "data_file_division": (
            left_record.annotations["data_file_division"]
            if ("data_file_division" in left_record.annotations)
            and (len(left_record.annotations["data_file_division"]) > 0)
            else (
                right_record.annotations["data_file_division"]
                if ("data_file_division" in right_record.annotations)
                and (len(right_record.annotations["data_file_division"]) > 0)
                else "UNK"
            )
        ),
        "date": datetime.now().strftime("%d-%b-%Y").upper(),
        "accessions": [f"{left_record.id}_{right_record.id}"],
    }

    for i in range(15, min(len(left_seq), len(right_seq))):
        if left_seq[-i:] == right_seq[:i]:
            break
    if i > 45:
        logger.warning(f"Overlapping region is too long: {i}")
    truncated_right_record = slice_seq_record_preserve_truncated(
        right_record, (i, len(right_seq))
    )
    linker_feature = SeqFeature(
        FeatureLocation(
            ExactPosition(len(left_seq) - i),
            ExactPosition(len(left_seq)),
            strand=None,
        ),
        type="misc_feature",
        qualifiers={
            "note": ["Gibson assembly linker"],
            "standard_name": ["Ligation"],
        },
    )
    left_record.features.append(linker_feature)
    merged_record = left_record + truncated_right_record
    merged_record.id = str(f"{left_record.id}_{right_record.id}")[:16]
    merged_record.annotations = annotations

    return merged_record
