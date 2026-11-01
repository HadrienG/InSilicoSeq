from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pysam
import pytest

from iss.generator import write_alignments


FR_COORDINATES = ("ref", 100, 105, 104, 109, 9, -9)


def make_read(read_id, sequence, qualities, cigar=None):
    read = SeqRecord(Seq(sequence), id=read_id)
    read.letter_annotations["phred_quality"] = qualities

    if cigar is not None:
        read.annotations["cigar"] = cigar

    return read


def write_and_read_alignments(
    tmp_path,
    coordinates,
    forward_cigar=None,
    reverse_cigar=None,
):
    sam_path = tmp_path / "test.sam"
    header = {
        "HD": {"VN": "1.6", "SO": "unsorted"},
        "SQ": [{"SN": "ref", "LN": 1000}],
    }

    forward = make_read(
        "pair/1",
        "ACGTA",
        [30, 31, 32, 33, 34],
        cigar=forward_cigar,
    )
    reverse = make_read(
        "pair/2",
        "TGCAT",
        [35, 36, 37, 38, 39],
        cigar=reverse_cigar,
    )

    with pysam.AlignmentFile(str(sam_path), "w", header=header) as sam:
        write_alignments(forward, reverse, coordinates, sam)

    with pysam.AlignmentFile(str(sam_path), "r") as sam:
        records = list(sam.fetch(until_eof=True))

    assert len(records) == 2
    read1 = next(record for record in records if record.is_read1)
    read2 = next(record for record in records if record.is_read2)
    return read1, read2


def test_write_alignments_sets_pair_fields(tmp_path):
    read1, read2 = write_and_read_alignments(tmp_path, FR_COORDINATES)

    assert read1.query_name == read2.query_name == "pair"
    assert read1.flag == 99
    assert read2.flag == 147

    assert read1.reference_start == 100
    assert read2.reference_start == 104
    assert read1.next_reference_start == 104
    assert read2.next_reference_start == 100


def test_write_alignments_orients_reverse_read_sequence_and_qualities(tmp_path):
    read1, read2 = write_and_read_alignments(tmp_path, FR_COORDINATES)

    assert read1.query_sequence == "ACGTA"
    assert read2.query_sequence == "ATGCA"
    assert read1.get_forward_sequence() == "ACGTA"
    assert read2.get_forward_sequence() == "TGCAT"

    assert read1.query_qualities.tolist() == [30, 31, 32, 33, 34]
    assert read2.query_qualities.tolist() == [39, 38, 37, 36, 35]


def test_write_alignments_uses_placeholder_match_cigar(tmp_path):
    read1, read2 = write_and_read_alignments(tmp_path, FR_COORDINATES)

    assert read1.cigarstring == "5M"
    assert read2.cigarstring == "5M"

    for record in (read1, read2):
        assert record.query_length == len(record.query_sequence)


def test_write_alignments_converts_insertion_annotation(tmp_path):
    read1, read2 = write_and_read_alignments(
        tmp_path,
        FR_COORDINATES,
        forward_cigar=[0, 0, 1, 0, 0],  # 2M1I2M
    )

    assert read1.cigarstring == "2M1I2M"
    assert read2.cigarstring == "5M"

    for record in (read1, read2):
        assert record.query_length == len(record.query_sequence)


@pytest.mark.parametrize(
    ("coordinates", "expected_read1_tlen", "expected_read2_tlen"),
    [
        pytest.param(
            ("ref", 100, 105, 104, 109, 9, -9),
            9,
            -9,
            id="read1-leftmost",
        ),
        pytest.param(
            ("ref", 200, 205, 190, 195, -15, 15),
            -15,
            15,
            id="fallback-read2-leftmost",
        ),
    ],
)
def test_write_alignments_assigns_tlen_by_leftmost_read(
    tmp_path, coordinates, expected_read1_tlen, expected_read2_tlen
):
    read1, read2 = write_and_read_alignments(tmp_path, coordinates)

    assert read1.template_length == expected_read1_tlen
    assert read2.template_length == expected_read2_tlen
