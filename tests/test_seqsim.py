"""
test_distance
=============

Tests for the `distance` module of the `seqsim` package.
"""

# Import Python standard libraries
import pytest

# Import the library being tested
import seqsim

test1_seq_a = "kitten"
test1_seq_b = "sitting"
test2_seq_a = (1, 2, 3)
test2_seq_b = [1, 2, 3]
test3_seq_a = (1, 2, 3, 4, 5)
test3_seq_b = (1, 2, 4, 3, 6, 7)
test4_seq_a = (1, 2, 3)
test4_seq_b = ["a", "b", "c", "d"]


def test_edit_distance():
    assert seqsim.distance.edit_distance(test1_seq_a, test1_seq_b) == pytest.approx(3.0)
    assert seqsim.distance.edit_distance(test2_seq_a, test2_seq_b) == pytest.approx(0.0)
    assert seqsim.distance.edit_distance(test3_seq_a, test3_seq_b) == pytest.approx(3.0)
    assert seqsim.distance.edit_distance(test4_seq_a, test4_seq_b) == pytest.approx(4.0)


def test_jaccard_distance():
    assert seqsim.distance.jaccard_distance(test1_seq_a, test1_seq_b) == pytest.approx(
        0.7
    )
    assert seqsim.distance.jaccard_distance(test2_seq_a, test2_seq_b) == pytest.approx(
        0.0
    )
    assert seqsim.distance.jaccard_distance(test3_seq_a, test3_seq_b) == pytest.approx(
        0.428571, rel=1e-5
    )
    assert seqsim.distance.jaccard_distance(test4_seq_a, test4_seq_b) == pytest.approx(
        1.0
    )


def test_subseq_jaccard_distance():
    assert seqsim.distance.subseq_jaccard_distance(
        test1_seq_a, test1_seq_b
    ) == pytest.approx(0.751556)
    assert seqsim.distance.subseq_jaccard_distance(
        test2_seq_a, test2_seq_b
    ) == pytest.approx(0.0)
    assert seqsim.distance.subseq_jaccard_distance(
        test3_seq_a, test3_seq_b
    ) == pytest.approx(0.787094)
    assert seqsim.distance.subseq_jaccard_distance(
        test4_seq_a, test4_seq_b
    ) == pytest.approx(1.0)


def test_mmcwpa_distance():
    assert seqsim.distance.mmcwpa_distance(test1_seq_a, test1_seq_b) == pytest.approx(
        0.538461, rel=1e-5
    )
    assert seqsim.distance.mmcwpa_distance(test2_seq_a, test2_seq_b) == pytest.approx(
        0.0
    )
    assert seqsim.distance.mmcwpa_distance(test3_seq_a, test3_seq_b) == pytest.approx(
        0.554638, rel=1e-5
    )
    assert seqsim.distance.mmcwpa_distance(test4_seq_a, test4_seq_b) == pytest.approx(
        1.0
    )
