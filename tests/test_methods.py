"""
test_distance
=============

Tests for the `distance` module of the `seqsim` package.
"""

# TODO: add tests with empty sequences

# Import Python standard libraries
import pytest

# Import the library being tested
from seqsim import methods

TEST1 = ["kitten", "sitting"]
TEST2 = [(1, 2, 3), [1, 2, 3]]
TEST3 = [(1, 2, 3, 4, 5), (1, 2, 4, 3, 6, 7)]
TEST4 = [(1, 2, 3), ["a", "b", "c", "d"]]

# TESTS FOR MEASURES OF SIMILARITY
# --------------------------------
#
# Measures of similarity don't offer all the properties of distance measures (i.e.,
# positivity, symmetry, identity-discerning, and triangle inequality), but may
# offer some of these.


@pytest.mark.parametrize(
    "seq_x,seq_y,expected,expected_norm,tol_norm",
    [
        TEST1 + [10.0, 0.3125, 1e-6],
        TEST2 + [6.0, 1.0, 0.0],
        TEST3 + [5.0, 0.238095, 1e-6],
        TEST4 + [0.0, 0.0, 0.0],
    ],
)
def test_birnbaum_simil(seq_x, seq_y, expected, expected_norm, tol_norm):
    # Test hard-coded expected value without normalization
    assert methods.birnbaum_simil(seq_x, seq_y) == expected

    # Test hard-coded expected value with normalization
    assert methods.birnbaum_simil(seq_x, seq_y, normal=True) == pytest.approx(
        expected_norm, abs=tol_norm
    )


@pytest.mark.parametrize(
    "seq_x,seq_y,expected,expected_norm,tol_norm",
    [
        TEST1 + [7.0, 0.25, 0.0],
        TEST2 + [6.0, 1.0, 0.0],
        TEST3 + [4.0, 0.190476, 1e-6],
        TEST4 + [0.0, 0.0, 0.0],
    ],
)
def test_fast_birnbaum_simil(seq_x, seq_y, expected, expected_norm, tol_norm):
    # Test hard-coded expected value without normalization
    assert methods.fast_birnbaum_simil(seq_x, seq_y) == expected

    # Test hard-coded expected value with normalization
    assert methods.fast_birnbaum_simil(seq_x, seq_y, normal=True) == pytest.approx(
        expected_norm, abs=tol_norm
    )


# TODO: add tests where fragile ends is different from levenshtein
# TODO: move fragile_end to similarity as it is not (necessarily) symmetric
# TODO: add a symmetric version with the mean?
@pytest.mark.parametrize(
    "seq_x,seq_y,expected,expected_norm,tol_norm",
    [
        TEST1 + [3.0, 0.5, 0.0],
        TEST2 + [0.0, 0.0, 0.0],
        TEST3 + [3.0, 0.6, 0.0],
        TEST4 + [4.0, 1.0, 0.0],
    ],
)
def test_fragile_ends_distance(seq_x, seq_y, expected, expected_norm, tol_norm):
    # Test hard-coded expected value without normalization
    assert methods.fragile_ends_simil(seq_x, seq_y) == pytest.approx(expected)

    # Test hard-coded expected value with normalization
    assert methods.fragile_ends_simil(seq_x, seq_y, normal=True) == pytest.approx(
        expected_norm, abs=tol_norm
    )

    # Test triangle-inequality
    seq_z = [element for element in seq_x] + [element for element in seq_y]
    dist_xz = methods.fragile_ends_simil(seq_x, seq_z)
    dist_xy = methods.fragile_ends_simil(seq_x, seq_y)
    dist_yz = methods.fragile_ends_simil(seq_y, seq_z)
    assert dist_xz <= (dist_xy + dist_yz)


# TODO: add tests where stemmatology is different from levenshtein
# TODO: move stemmatological to similarity as it is not (necessarily) symmetric
# TODO: add a symmetric version with the mean?
@pytest.mark.parametrize(
    "seq_x,seq_y,expected,expected_norm,tol_norm",
    [
        TEST1 + [3.0, 0.428571, 1e-6],
        TEST2 + [0.0, 0.0, 0.0],
        TEST3 + [3.0, 0.5, 0.0],
        TEST4 + [4.0, 1.0, 0.0],
    ],
)
def test_stemmatology_distance(seq_x, seq_y, expected, expected_norm, tol_norm):
    # Test hard-coded expected value without normalization
    assert methods.stemmatological_simil(seq_x, seq_y) == expected

    # Test hard-coded expected value without normalization
    assert methods.stemmatological_simil(seq_x, seq_y, normal=True) == pytest.approx(
        expected_norm, abs=tol_norm
    )

    # Test triangle-inequality
    seq_z = [element for element in seq_x] + [element for element in seq_y]
    dist_xz = methods.stemmatological_simil(seq_x, seq_z)
    dist_xy = methods.stemmatological_simil(seq_x, seq_y)
    dist_yz = methods.stemmatological_simil(seq_y, seq_z)
    assert dist_xz <= (dist_xy + dist_yz)


# TODO: add tests where stemmatology (norm/2030) is different from levenshtein (norm)
# TODO: move norm_stemmatological to similarity as it is not (necessarily) symmetric
# TODO: add a symmetric version with the mean?
# TODO: add non-normal test (if we keep this)
@pytest.mark.parametrize(
    "seq_x,seq_y,expected,tol",
    [
        TEST1 + [0.428571, 1e-6],
        TEST2 + [0.0, 0.0],
        TEST3 + [0.5, 0.0],
        TEST4 + [1.0, 0.0],
    ],
)
def test_stemmatology_2030_distance(seq_x, seq_y, expected, tol):
    # Test hard-coded expected value
    assert methods.stemmatological_simil(
        seq_x, seq_y, 20.0, 30.0, normal=True
    ) == pytest.approx(expected, abs=tol)

    # Test triangle-inequality
    seq_z = [element for element in seq_x] + [element for element in seq_y]
    dist_xz = methods.stemmatological_simil(seq_x, seq_z, 20.0, 30.0, normal=True)
    dist_xy = methods.stemmatological_simil(seq_x, seq_y, 20.0, 30.0, normal=True)
    dist_yz = methods.stemmatological_simil(seq_y, seq_z, 20.0, 30.0, normal=True)
    assert dist_xz <= (dist_xy + dist_yz)


# TESTS FOR MEASURES OF DISTANCE
# ------------------------------
#
# Measures of distance offer all the properties of distance measures (i.e.,


@pytest.mark.parametrize(
    "seq_x,seq_y,expected,tol",
    [
        TEST1 + [0.7, 0.0],
        TEST2 + [0.0, 0.0],
        TEST3 + [0.428571, 1e-6],
        TEST4 + [1.0, 0.0],
    ],
)
def test_jaccard_distance(seq_x, seq_y, expected, tol):
    # Test hard-coded expected value
    assert methods.jaccard_dist(seq_x, seq_y) == pytest.approx(expected, abs=tol)

    # Test symmetry
    assert methods.jaccard_dist(seq_x, seq_y) == methods.jaccard_dist(seq_y, seq_x)

    # Test triangle-inequality
    seq_z = [element for element in seq_x] + [element for element in seq_y]
    dist_xz = methods.jaccard_dist(seq_x, seq_z)
    dist_xy = methods.jaccard_dist(seq_x, seq_y)
    dist_yz = methods.jaccard_dist(seq_y, seq_z)
    assert dist_xz <= (dist_xy + dist_yz)


@pytest.mark.parametrize(
    "seq_x,seq_y,expected,tol",
    [
        TEST1 + [0.751556, 1e-6],
        TEST2 + [0.0, 0.0],
        TEST3 + [0.787094, 1e-6],
        TEST4 + [1.0, 0.0],
    ],
)
def test_subseq_jaccard_distance(seq_x, seq_y, expected, tol):
    # Test hard-coded expected value
    assert methods.subseq_jaccard_dist(seq_x, seq_y) == pytest.approx(expected, abs=tol)

    # Test symmetry
    assert methods.subseq_jaccard_dist(seq_x, seq_y) == methods.subseq_jaccard_dist(
        seq_y, seq_x
    )

    # Test triangle-inequality
    seq_z = [element for element in seq_x] + [element for element in seq_y]
    dist_xz = methods.subseq_jaccard_dist(seq_x, seq_z)
    dist_xy = methods.subseq_jaccard_dist(seq_x, seq_y)
    dist_yz = methods.subseq_jaccard_dist(seq_y, seq_z)
    assert dist_xz <= (dist_xy + dist_yz)


@pytest.mark.parametrize(
    "seq_x,seq_y,expected,expected_norm,tol_norm",
    [
        TEST1 + [3.0, 0.428571, 1e-6],
        TEST2 + [0.0, 0.0, 0.0],
        TEST3 + [3.0, 0.5, 0.0],
        TEST4 + [4.0, 1.0, 0.0],
    ],
)
def test_levenshtein_distance(seq_x, seq_y, expected, expected_norm, tol_norm):
    # Test hard-coded expected value without normalization
    assert methods.levenshtein_dist(seq_x, seq_y) == expected

    # Test hard-coded expected value with normalization
    assert methods.levenshtein_dist(seq_x, seq_y, normal=True) == pytest.approx(
        expected_norm, abs=tol_norm
    )

    # Test symmetry
    assert methods.levenshtein_dist(seq_x, seq_y) == methods.levenshtein_dist(
        seq_y, seq_x
    )

    # Test triangle-inequality
    seq_z = [element for element in seq_x] + [element for element in seq_y]
    dist_xz = methods.levenshtein_dist(seq_x, seq_z)
    dist_xy = methods.levenshtein_dist(seq_x, seq_y)
    dist_yz = methods.levenshtein_dist(seq_y, seq_z)
    assert dist_xz <= (dist_xy + dist_yz)


# TODO: add tests where levdamerau is different from levenshtein
@pytest.mark.parametrize(
    "seq_x,seq_y,expected,expected_norm,tol_norm",
    [
        TEST1 + [3.0, 0.428571, 1e-6],
        TEST2 + [0.0, 0.0, 0.0],
        TEST3 + [3.0, 0.5, 0.0],
        TEST4 + [4.0, 1.0, 0.0],
    ],
)
def test_levdamerau_distance(seq_x, seq_y, expected, expected_norm, tol_norm):
    # Test hard-coded expected value without normalization
    assert methods.levdamerau_dist(seq_x, seq_y) == expected

    # Test hard-coded expected value with normalization
    assert methods.levdamerau_dist(seq_x, seq_y, normal=True) == pytest.approx(
        expected_norm, abs=tol_norm
    )

    # Test symmetry
    assert methods.levdamerau_dist(seq_x, seq_y) == methods.levdamerau_dist(
        seq_y, seq_x
    )

    # Test triangle-inequality
    seq_z = [element for element in seq_x] + [element for element in seq_y]
    dist_xz = methods.levdamerau_dist(seq_x, seq_z)
    dist_xy = methods.levdamerau_dist(seq_x, seq_y)
    dist_yz = methods.levdamerau_dist(seq_y, seq_z)
    assert dist_xz <= (dist_xy + dist_yz)


# TODO: add tests where bulk delete is different from levenshtein
@pytest.mark.parametrize(
    "seq_x,seq_y,expected,expected_norm,tol_norm",
    [
        TEST1 + [3.0, 0.428571, 1e-6],
        TEST2 + [0.0, 0.0, 0.0],
        TEST3 + [3.0, 0.5, 0.0],
        TEST4 + [4.0, 1.0, 0.0],
    ],
)
def test_bulk_delete_distance(seq_x, seq_y, expected, expected_norm, tol_norm):
    # Test hard-coded expected value without normalization
    assert methods.bulk_delete_dist(seq_x, seq_y) == expected

    # Test hard-coded expected value without normalization
    assert methods.bulk_delete_dist(seq_x, seq_y, normal=True) == pytest.approx(
        expected_norm, abs=tol_norm
    )

    # Test symmetry
    assert methods.bulk_delete_dist(seq_x, seq_y) == methods.bulk_delete_dist(
        seq_y, seq_x
    )

    # Test triangle-inequality
    seq_z = [element for element in seq_x] + [element for element in seq_y]
    dist_xz = methods.bulk_delete_dist(seq_x, seq_z)
    dist_xy = methods.bulk_delete_dist(seq_x, seq_y)
    dist_yz = methods.bulk_delete_dist(seq_y, seq_z)
    assert dist_xz <= (dist_xy + dist_yz)


@pytest.mark.parametrize(
    "seq_x,seq_y,expected,tol",
    [
        TEST1 + [0.253968, 1e-6],
        TEST2 + [0.0, 0.0],
        TEST3 + [0.261111, 1e-6],
        TEST4 + [1.0, 0.0],
    ],
)
def test_jaro_distance(seq_x, seq_y, expected, tol):
    # Test hard-coded expected value
    assert methods.jaro_dist(seq_x, seq_y) == pytest.approx(expected, abs=tol)

    # Test symmetry
    assert methods.jaro_dist(seq_x, seq_y) == methods.jaro_dist(seq_y, seq_x)

    # Test triangle-inequality
    seq_z = [element for element in seq_x] + [element for element in seq_y]
    dist_xz = methods.jaro_dist(seq_x, seq_z)
    dist_xy = methods.jaro_dist(seq_x, seq_y)
    dist_yz = methods.jaro_dist(seq_y, seq_z)
    assert dist_xz <= (dist_xy + dist_yz)


@pytest.mark.parametrize(
    "seq_x,seq_y,expected,tol",
    [
        TEST1 + [0.253968, 1e-6],
        TEST2 + [0.0, 0.0],
        TEST3 + [0.208888, 1e-6],
        TEST4 + [1.0, 0.0],
    ],
)
def test_jarowinkler_distance(seq_x, seq_y, expected, tol):
    # Test hard-coded expected value
    assert methods.jaro_winkler_dist(seq_x, seq_y) == pytest.approx(expected, abs=tol)

    # Test symmetry
    assert methods.jaro_winkler_dist(seq_x, seq_y) == methods.jaro_winkler_dist(
        seq_y, seq_x
    )

    # Test triangle-inequality
    seq_z = [element for element in seq_x] + [element for element in seq_y]
    dist_xz = methods.jaro_winkler_dist(seq_x, seq_z)
    dist_xy = methods.jaro_winkler_dist(seq_x, seq_y)
    dist_yz = methods.jaro_winkler_dist(seq_y, seq_z)
    assert dist_xz <= (dist_xy + dist_yz)


@pytest.mark.parametrize(
    "seq_x,seq_y,expected,tol",
    [
        TEST1 + [0.538461, 1e-6],
        TEST2 + [0.0, 0.0],
        TEST3 + [0.554638, 1e-6],
        TEST4 + [1.0, 0.0],
    ],
)
def test_mmcwpa_distance(seq_x, seq_y, expected, tol):
    # Test hard-coded expected value
    assert methods.mmcwpa_dist(seq_x, seq_y) == pytest.approx(expected, abs=tol)

    # Test symmetry
    assert methods.mmcwpa_dist(seq_x, seq_y) == methods.mmcwpa_dist(seq_y, seq_x)

    # Test triangle-inequality
    seq_z = [element for element in seq_x] + [element for element in seq_y]
    dist_xz = methods.mmcwpa_dist(seq_x, seq_z)
    dist_xy = methods.mmcwpa_dist(seq_x, seq_y)
    dist_yz = methods.mmcwpa_dist(seq_y, seq_z)
    assert dist_xz <= (dist_xy + dist_yz)


@pytest.mark.parametrize(
    "seq_x,seq_y,expected,tol",
    [
        TEST1 + [0.565217, 1e-6],
        TEST2 + [0.0, 0.0],
        TEST3 + [0.666666, 1e-6],
        TEST4 + [1.0, 0.0],
    ],
)
def test_birnbaum_distance(seq_x, seq_y, expected, tol):
    # Test hard-coded expected value
    assert methods.birnbaum_dist(seq_x, seq_y) == pytest.approx(expected, abs=tol)

    # Test symmetry
    assert methods.birnbaum_dist(seq_x, seq_y) == methods.birnbaum_dist(seq_y, seq_x)

    # Test triangle-inequality
    seq_z = [element for element in seq_x] + [element for element in seq_y]
    dist_xz = methods.birnbaum_dist(seq_x, seq_z)
    dist_xy = methods.birnbaum_dist(seq_x, seq_y)
    dist_yz = methods.birnbaum_dist(seq_y, seq_z)
    assert dist_xz <= (dist_xy + dist_yz)


@pytest.mark.parametrize(
    "seq_x,seq_y,expected,tol",
    [
        TEST1 + [0.666666, 1e-6],
        TEST2 + [0.0, 0.0],
        TEST3 + [0.733333, 1e-6],
        TEST4 + [1.0, 0.0],
    ],
)
def test_fast_birnbaum_distance(seq_x, seq_y, expected, tol):
    # Test hard-coded expected value
    assert methods.fast_birnbaum_dist(seq_x, seq_y) == pytest.approx(expected, abs=tol)

    # Test symmetry
    assert methods.fast_birnbaum_dist(seq_x, seq_y) == methods.fast_birnbaum_dist(
        seq_y, seq_x
    )

    # Test triangle-inequality
    seq_z = [element for element in seq_x] + [element for element in seq_y]
    dist_xz = methods.fast_birnbaum_dist(seq_x, seq_z)
    dist_xy = methods.fast_birnbaum_dist(seq_x, seq_y)
    dist_yz = methods.fast_birnbaum_dist(seq_y, seq_z)
    assert dist_xz <= (dist_xy + dist_yz)
