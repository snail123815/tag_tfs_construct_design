
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import MeltingTemp as mt

def find_local_best_primer_bind(
    primer_extended: Seq | SeqRecord,
    target_primer_size: int = 20,
    wander_range: tuple[int] = (-1, 6),
    tm_range: tuple[float] = (58, 68),
    panalty_offset: int = 1,
    panalty_3p_repeat: int = 5,
    panalty_tm_not_in_range: int = 10,
):
    # Check differences around the 20th nucleotide for repeats
    if isinstance(primer_extended, SeqRecord):
        primer_extended = primer_extended.seq
    target_primer = primer_extended[:target_primer_size]
    assert (
        len(primer_extended) >= target_primer_size + wander_range[1]
    ), "Provided sequence (primer_extended) must be enough for wander_range"

    putatives = []
    assert wander_range[0] > -6, "Wander must start after -6"
    assert wander_range[1] <= 14, "Wander must end before 14"
    wander_offsets = []
    for offset in [0, 1, 2, -1, 3, 4, -2, 5, 6, -3, 7, 8] + [
        -4,
        9,
        10,
        11,
        -5,
        12,
        13,
        14,
        -6,
    ]:
        if offset in range(wander_range[0], wander_range[1]):
            wander_offsets.append(offset)

    for offset_panalty, offset in enumerate(wander_offsets):
        score = -offset_panalty * panalty_offset
        tp = primer_extended[: target_primer_size + offset]
        if tp[-2] == tp[-1]:
            score -= panalty_3p_repeat
        if mt.Tm_NN(tp) < tm_range[0] or mt.Tm_NN(tp) > tm_range[1]:
            score -= panalty_tm_not_in_range
        putatives.append((tp, score))
    target_primer, score = max(putatives, key=lambda x: x[1])

    return target_primer, score
