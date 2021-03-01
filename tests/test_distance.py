"""
test_distance
=============

Tests for the `distance` module of the `seqsim` package.
"""

# TODO: add tests confirming the distance properties: positivity, symmetry, id, tri-ineq

# Import Python standard libraries
import pytest

# Import the library being tested
from seqsim import distance

test1_a = "kitten"
test1_b = "sitting"
test2_a = (1, 2, 3)
test2_b = [1, 2, 3]
test3_a = (1, 2, 3, 4, 5)
test3_b = (1, 2, 4, 3, 6, 7)
test4_a = (1, 2, 3)
test4_b = ["a", "b", "c", "d"]


def test_birnbaum_distance():
    assert distance.fast_birnbaum(test1_a, test1_b) == pytest.approx(0.666666, abs=1e-6)
    assert distance.fast_birnbaum(test2_a, test2_b) == pytest.approx(0.0)
    assert distance.fast_birnbaum(test3_a, test3_b) == pytest.approx(0.733333, abs=1e-6)
    assert distance.fast_birnbaum(test4_a, test4_b) == pytest.approx(1.0, abs=1e-6)


def test_levenshtein_distance():
    assert distance.levenshtein(test1_a, test1_b) == pytest.approx(3.0)
    assert distance.levenshtein(test2_a, test2_b) == pytest.approx(0.0)
    assert distance.levenshtein(test3_a, test3_b) == pytest.approx(3.0)
    assert distance.levenshtein(test4_a, test4_b) == pytest.approx(4.0)


def test_norm_levenshtein_distance():
    assert distance.norm_levenshtein(test1_a, test1_b) == pytest.approx(
        0.428571, abs=1e-6
    )
    assert distance.norm_levenshtein(test2_a, test2_b) == pytest.approx(0.0)
    assert distance.norm_levenshtein(test3_a, test3_b) == pytest.approx(0.5)
    assert distance.norm_levenshtein(test4_a, test4_b) == pytest.approx(1.0)


# TODO: add tests where levdamerau is different from levenshtein
def test_levdamerau_distance():
    assert distance.levdamerau(test1_a, test1_b) == pytest.approx(3.0)
    assert distance.levdamerau(test2_a, test2_b) == pytest.approx(0.0)
    assert distance.levdamerau(test3_a, test3_b) == pytest.approx(3.0)
    assert distance.levdamerau(test4_a, test4_b) == pytest.approx(4.0)


# TODO: add tests where fragile ends is different from levenshtein
def test_fragile_ends_distance():
    assert distance.fragile_ends(test1_a, test1_b) == pytest.approx(3.0)
    assert distance.fragile_ends(test2_a, test2_b) == pytest.approx(0.0)
    assert distance.fragile_ends(test3_a, test3_b) == pytest.approx(3.0)
    assert distance.fragile_ends(test4_a, test4_b) == pytest.approx(4.0)


# TODO: add tests where bulk delete is different from levenshtein
def test_bulk_delete_distance():
    assert distance.bulk_delete(test1_a, test1_b) == pytest.approx(3.0)
    assert distance.bulk_delete(test2_a, test2_b) == pytest.approx(0.0)
    assert distance.bulk_delete(test3_a, test3_b) == pytest.approx(3.0)
    assert distance.bulk_delete(test4_a, test4_b) == pytest.approx(4.0)


# TODO: add tests where stemmatology is different from levenshtein
def test_stemmatology_distance():
    assert distance.stemmatological(test1_a, test1_b) == pytest.approx(3.0)
    assert distance.stemmatological(test2_a, test2_b) == pytest.approx(0.0)
    assert distance.stemmatological(test3_a, test3_b) == pytest.approx(3.0)
    assert distance.stemmatological(test4_a, test4_b) == pytest.approx(4.0)


# TODO: add tests where stemmatology (norm) is different from levenshtein (norm)
def test_norm_stemmatology_distance():
    assert distance.norm_stemmatological(test1_a, test1_b) == pytest.approx(
        0.428571, abs=1e-6
    )
    assert distance.norm_stemmatological(test2_a, test2_b) == pytest.approx(0.0)
    assert distance.norm_stemmatological(test3_a, test3_b) == pytest.approx(0.5)
    assert distance.norm_stemmatological(test4_a, test4_b) == pytest.approx(1.0)


# TODO: add tests where stemmatology (norm/2030) is different from levenshtein (norm)
def test_norm_stemmatology_2030_distance():
    assert distance.norm_stemmatological_2030(test1_a, test1_b) == pytest.approx(
        0.428571, abs=1e-6
    )
    assert distance.norm_stemmatological_2030(test2_a, test2_b) == pytest.approx(0.0)
    assert distance.norm_stemmatological_2030(test3_a, test3_b) == pytest.approx(0.5)
    assert distance.norm_stemmatological_2030(test4_a, test4_b) == pytest.approx(1.0)


def test_jaccard_distance():
    assert distance.jaccard(test1_a, test1_b) == pytest.approx(0.7)
    assert distance.jaccard(test2_a, test2_b) == pytest.approx(0.0)
    assert distance.jaccard(test3_a, test3_b) == pytest.approx(0.428571, rel=1e-5)
    assert distance.jaccard(test4_a, test4_b) == pytest.approx(1.0)


def test_subseq_jaccard_distance():
    assert distance.subseq_jaccard(test1_a, test1_b) == pytest.approx(0.751556)
    assert distance.subseq_jaccard(test2_a, test2_b) == pytest.approx(0.0)
    assert distance.subseq_jaccard(test3_a, test3_b) == pytest.approx(0.787094)
    assert distance.subseq_jaccard(test4_a, test4_b) == pytest.approx(1.0)


def test_mmcwpa_distance():
    assert distance.mmcwpa(test1_a, test1_b) == pytest.approx(0.538461, rel=1e-5)
    assert distance.mmcwpa(test2_a, test2_b) == pytest.approx(0.0)
    assert distance.mmcwpa(test3_a, test3_b) == pytest.approx(0.554638, rel=1e-5)
    assert distance.mmcwpa(test4_a, test4_b) == pytest.approx(1.0)
