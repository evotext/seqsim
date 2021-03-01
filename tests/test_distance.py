"""
test_distance
=============

Tests for the `distance` module of the `seqsim` package.
"""

# TODO: add tests confirming the distance properties: positivity, symmetry, id, tri-ineq
# TODO: add tests with empty sequences

# Import Python standard libraries
import pytest

# Import the library being tested
from seqsim import distance

test1 = ["kitten", "sitting"]
test2 = [(1, 2, 3), [1, 2, 3]]
test3 = [(1, 2, 3, 4, 5), (1, 2, 4, 3, 6, 7)]
test4 = [(1, 2, 3), ["a", "b", "c", "d"]]

test1_a = "kitten"
test1_b = "sitting"
test2_a = (1, 2, 3)
test2_b = [1, 2, 3]
test3_a = (1, 2, 3, 4, 5)
test3_b = (1, 2, 4, 3, 6, 7)
test4_a = (1, 2, 3)
test4_b = ["a", "b", "c", "d"]


@pytest.mark.parametrize(
    "seq_x,seq_y,expected,tol",
    [
        test1 + [0.666666, 1e-6],
        test2 + [0.0, 0.0],
        test3 + [0.733333, 1e-6],
        test4 + [1.0, 0.0],
    ],
)
def test_fast_birnbaum_distance(seq_x, seq_y, expected, tol):
    # Test hard-coded expected value
    assert distance.fast_birnbaum(seq_x, seq_y) == pytest.approx(expected, abs=tol)

    # Test symmetry
    assert distance.fast_birnbaum(seq_x, seq_y) == distance.fast_birnbaum(seq_y, seq_x)

    # Test triangle-inequality
    seq_z = [element for element in seq_x] + [element for element in seq_y]
    dist_xz = distance.fast_birnbaum(seq_x, seq_z)
    dist_xy = distance.fast_birnbaum(seq_x, seq_y)
    dist_yz = distance.fast_birnbaum(seq_y, seq_z)
    assert dist_xz <= (dist_xy + dist_yz)


@pytest.mark.parametrize(
    "seq_x,seq_y,expected,tol",
    [
        test1 + [3.0, 0.0],
        test2 + [0.0, 0.0],
        test3 + [3.0, 0.0],
        test4 + [4.0, 0.0],
    ],
)
def test_levenshtein_distance(seq_x, seq_y, expected, tol):
    # Test hard-coded expected value
    assert distance.levenshtein(seq_x, seq_y) == pytest.approx(expected, abs=tol)

    # Test symmetry
    assert distance.levenshtein(seq_x, seq_y) == distance.levenshtein(seq_y, seq_x)

    # Test triangle-inequality
    seq_z = [element for element in seq_x] + [element for element in seq_y]
    dist_xz = distance.levenshtein(seq_x, seq_z)
    dist_xy = distance.levenshtein(seq_x, seq_y)
    dist_yz = distance.levenshtein(seq_y, seq_z)
    assert dist_xz <= (dist_xy + dist_yz)


@pytest.mark.parametrize(
    "seq_x,seq_y,expected,tol",
    [
        test1 + [0.428571, 1e-6],
        test2 + [0.0, 0.0],
        test3 + [0.5, 0.0],
        test4 + [1.0, 0.0],
    ],
)
def test_norm_levenshtein_distance(seq_x, seq_y, expected, tol):
    # Test hard-coded expected value
    assert distance.norm_levenshtein(seq_x, seq_y) == pytest.approx(expected, abs=tol)

    # Test symmetry
    assert distance.norm_levenshtein(seq_x, seq_y) == distance.norm_levenshtein(
        seq_y, seq_x
    )

    # Test triangle-inequality
    seq_z = [element for element in seq_x] + [element for element in seq_y]
    dist_xz = distance.norm_levenshtein(seq_x, seq_z)
    dist_xy = distance.norm_levenshtein(seq_x, seq_y)
    dist_yz = distance.norm_levenshtein(seq_y, seq_z)
    assert dist_xz <= (dist_xy + dist_yz)


# TODO: add tests where levdamerau is different from levenshtein
@pytest.mark.parametrize(
    "seq_x,seq_y,expected,tol",
    [
        test1 + [3.0, 0.0],
        test2 + [0.0, 0.0],
        test3 + [3.0, 0.0],
        test4 + [4.0, 0.0],
    ],
)
def test_levdamerau_distance(seq_x, seq_y, expected, tol):
    # Test hard-coded expected value
    assert distance.levdamerau(seq_x, seq_y) == pytest.approx(expected, abs=tol)

    # Test symmetry
    assert distance.levdamerau(seq_x, seq_y) == distance.levdamerau(seq_y, seq_x)

    # Test triangle-inequality
    seq_z = [element for element in seq_x] + [element for element in seq_y]
    dist_xz = distance.levdamerau(seq_x, seq_z)
    dist_xy = distance.levdamerau(seq_x, seq_y)
    dist_yz = distance.levdamerau(seq_y, seq_z)
    assert dist_xz <= (dist_xy + dist_yz)


# TODO: add tests where fragile ends is different from levenshtein
# TODO: move fragile_end to similarity as it is not (necessarily) symmetric
# TODO: add a symmetric version with the mean?
@pytest.mark.parametrize(
    "seq_x,seq_y,expected,tol",
    [
        test1 + [3.0, 0.0],
        test2 + [0.0, 0.0],
        test3 + [3.0, 0.0],
        test4 + [4.0, 0.0],
    ],
)
def test_fragile_ends_distance(seq_x, seq_y, expected, tol):
    # Test hard-coded expected value
    assert distance.fragile_ends(seq_x, seq_y) == pytest.approx(expected, abs=tol)

    # Test symmetry
    #    assert distance.fragile_ends(seq_x, seq_y) == distance.fragile_ends(seq_y, seq_x)

    # Test triangle-inequality
    seq_z = [element for element in seq_x] + [element for element in seq_y]
    dist_xz = distance.fragile_ends(seq_x, seq_z)
    dist_xy = distance.fragile_ends(seq_x, seq_y)
    dist_yz = distance.fragile_ends(seq_y, seq_z)
    assert dist_xz <= (dist_xy + dist_yz)


# TODO: add tests where bulk delete is different from levenshtein
@pytest.mark.parametrize(
    "seq_x,seq_y,expected,tol",
    [
        test1 + [3.0, 0.0],
        test2 + [0.0, 0.0],
        test3 + [3.0, 0.0],
        test4 + [4.0, 0.0],
    ],
)
def test_bulk_delete_distance(seq_x, seq_y, expected, tol):
    # Test hard-coded expected value
    assert distance.bulk_delete(seq_x, seq_y) == pytest.approx(expected, abs=tol)

    # Test symmetry
    assert distance.bulk_delete(seq_x, seq_y) == distance.bulk_delete(seq_y, seq_x)

    # Test triangle-inequality
    seq_z = [element for element in seq_x] + [element for element in seq_y]
    dist_xz = distance.bulk_delete(seq_x, seq_z)
    dist_xy = distance.bulk_delete(seq_x, seq_y)
    dist_yz = distance.bulk_delete(seq_y, seq_z)
    assert dist_xz <= (dist_xy + dist_yz)


# TODO: add tests where stemmatology is different from levenshtein
# TODO: move stemmatological to similarity as it is not (necessarily) symmetric
# TODO: add a symmetric version with the mean?
@pytest.mark.parametrize(
    "seq_x,seq_y,expected,tol",
    [
        test1 + [3.0, 0.0],
        test2 + [0.0, 0.0],
        test3 + [3.0, 0.0],
        test4 + [4.0, 0.0],
    ],
)
def test_stemmatology_distance(seq_x, seq_y, expected, tol):
    # Test hard-coded expected value
    assert distance.stemmatological(seq_x, seq_y) == pytest.approx(expected, abs=tol)

    # Test symmetry
    #    assert distance.stemmatological(seq_x, seq_y) == distance.stemmatological(seq_y, seq_x)

    # Test triangle-inequality
    seq_z = [element for element in seq_x] + [element for element in seq_y]
    dist_xz = distance.stemmatological(seq_x, seq_z)
    dist_xy = distance.stemmatological(seq_x, seq_y)
    dist_yz = distance.stemmatological(seq_y, seq_z)
    assert dist_xz <= (dist_xy + dist_yz)


# TODO: add tests where stemmatology (norm) is different from levenshtein (norm)
# TODO: move norm_stemmatological to similarity as it is not (necessarily) symmetric
# TODO: add a symmetric version with the mean?
@pytest.mark.parametrize(
    "seq_x,seq_y,expected,tol",
    [
        test1 + [0.428571, 1e-6],
        test2 + [0.0, 0.0],
        test3 + [0.5, 0.0],
        test4 + [1.0, 0.0],
    ],
)
def test_norm_stemmatology_distance(seq_x, seq_y, expected, tol):
    # Test hard-coded expected value
    assert distance.norm_stemmatological(seq_x, seq_y) == pytest.approx(
        expected, abs=tol
    )

    # Test symmetry
    #    assert distance.norm_stemmatological(seq_x, seq_y) == distance.norm_stemmatological(seq_y, seq_x)

    # Test triangle-inequality
    seq_z = [element for element in seq_x] + [element for element in seq_y]
    dist_xz = distance.norm_stemmatological(seq_x, seq_z)
    dist_xy = distance.norm_stemmatological(seq_x, seq_y)
    dist_yz = distance.norm_stemmatological(seq_y, seq_z)
    assert dist_xz <= (dist_xy + dist_yz)


# TODO: add tests where stemmatology (norm/2030) is different from levenshtein (norm)
# TODO: move norm_stemmatological to similarity as it is not (necessarily) symmetric
# TODO: add a symmetric version with the mean?
@pytest.mark.parametrize(
    "seq_x,seq_y,expected,tol",
    [
        test1 + [0.428571, 1e-6],
        test2 + [0.0, 0.0],
        test3 + [0.5, 0.0],
        test4 + [1.0, 0.0],
    ],
)
def test_norm_stemmatology_2030_distance(seq_x, seq_y, expected, tol):
    # Test hard-coded expected value
    assert distance.norm_stemmatological_2030(seq_x, seq_y) == pytest.approx(
        expected, abs=tol
    )

    # Test symmetry
    #    assert distance.norm_stemmatological_2030(seq_x, seq_y) == distance.norm_stemmatological_2030(seq_y, seq_x)

    # Test triangle-inequality
    seq_z = [element for element in seq_x] + [element for element in seq_y]
    dist_xz = distance.norm_stemmatological_2030(seq_x, seq_z)
    dist_xy = distance.norm_stemmatological_2030(seq_x, seq_y)
    dist_yz = distance.norm_stemmatological_2030(seq_y, seq_z)
    assert dist_xz <= (dist_xy + dist_yz)


@pytest.mark.parametrize(
    "seq_x,seq_y,expected,tol",
    [
        test1 + [0.7, 0.0],
        test2 + [0.0, 0.0],
        test3 + [0.428571, 1e-6],
        test4 + [1.0, 0.0],
    ],
)
def test_jaccard_distance(seq_x, seq_y, expected, tol):
    # Test hard-coded expected value
    assert distance.jaccard(seq_x, seq_y) == pytest.approx(expected, abs=tol)

    # Test symmetry
    assert distance.jaccard(seq_x, seq_y) == distance.jaccard(seq_y, seq_x)

    # Test triangle-inequality
    seq_z = [element for element in seq_x] + [element for element in seq_y]
    dist_xz = distance.jaccard(seq_x, seq_z)
    dist_xy = distance.jaccard(seq_x, seq_y)
    dist_yz = distance.jaccard(seq_y, seq_z)
    assert dist_xz <= (dist_xy + dist_yz)


@pytest.mark.parametrize(
    "seq_x,seq_y,expected,tol",
    [
        test1 + [0.751556, 1e-6],
        test2 + [0.0, 0.0],
        test3 + [0.787094, 1e-6],
        test4 + [1.0, 0.0],
    ],
)
def test_subseq_jaccard_distance(seq_x, seq_y, expected, tol):
    # Test hard-coded expected value
    assert distance.subseq_jaccard(seq_x, seq_y) == pytest.approx(expected, abs=tol)

    # Test symmetry
    assert distance.subseq_jaccard(seq_x, seq_y) == distance.subseq_jaccard(
        seq_y, seq_x
    )

    # Test triangle-inequality
    seq_z = [element for element in seq_x] + [element for element in seq_y]
    dist_xz = distance.subseq_jaccard(seq_x, seq_z)
    dist_xy = distance.subseq_jaccard(seq_x, seq_y)
    dist_yz = distance.subseq_jaccard(seq_y, seq_z)
    assert dist_xz <= (dist_xy + dist_yz)


@pytest.mark.parametrize(
    "seq_x,seq_y,expected,tol",
    [
        test1 + [0.538461, 1e-6],
        test2 + [0.0, 0.0],
        test3 + [0.554638, 1e-6],
        test4 + [1.0, 0.0],
    ],
)
def test_mmcwpa_distance(seq_x, seq_y, expected, tol):
    # Test hard-coded expected value
    assert distance.mmcwpa(seq_x, seq_y) == pytest.approx(expected, abs=tol)

    # Test symmetry
    assert distance.mmcwpa(seq_x, seq_y) == distance.mmcwpa(seq_y, seq_x)

    # Test triangle-inequality
    seq_z = [element for element in seq_x] + [element for element in seq_y]
    dist_xz = distance.mmcwpa(seq_x, seq_z)
    dist_xy = distance.mmcwpa(seq_x, seq_y)
    dist_yz = distance.mmcwpa(seq_y, seq_z)
    assert dist_xz <= (dist_xy + dist_yz)
