"""
test_edit
=========

Tests for the `edit` module of the `seqsim` package.
"""

# TODO: add tests with empty sequences

# Import Python standard libraries
import pytest

# Import the library being tested
from seqsim import edit

# TESTS FOR MEASURES OF SIMILARITY
# --------------------------------
#
# Measures of similarity don't offer all the properties of distance measures (i.e.,
# positivity, symmetry, identity-discerning, and triangle inequality), but may
# offer some of these.


@pytest.mark.parametrize(
    "seq_x,seq_y,expected,expected_norm,tol_norm",
    [
        ["kitten", "sitting"] + [10.0, 0.3125, 1e-6],
        [(1, 2, 3), [1, 2, 3]] + [6.0, 1.0, 0.0],
        [(1, 2, 3, 4, 5), (1, 2, 4, 3, 6, 7)] + [5.0, 0.238095, 1e-6],
        [(1, 2, 3), ["a", "b", "c", "d"]] + [0.0, 0.0, 0.0],
    ],
)
def test_birnbaum_simil(seq_x, seq_y, expected, expected_norm, tol_norm):
    # Test hard-coded expected value without normalization
    assert edit.birnbaum_simil(seq_x, seq_y) == expected

    # Test hard-coded expected value with normalization
    assert edit.birnbaum_simil(seq_x, seq_y, normal=True) == pytest.approx(
        expected_norm, abs=tol_norm
    )


@pytest.mark.parametrize(
    "seq_x,seq_y,expected,expected_norm,tol_norm",
    [
        ["kitten", "sitting"] + [7.0, 0.25, 0.0],
        [(1, 2, 3), [1, 2, 3]] + [6.0, 1.0, 0.0],
        [(1, 2, 3, 4, 5), (1, 2, 4, 3, 6, 7)] + [4.0, 0.190476, 1e-6],
        [(1, 2, 3), ["a", "b", "c", "d"]] + [0.0, 0.0, 0.0],
    ],
)
def test_fast_birnbaum_simil(seq_x, seq_y, expected, expected_norm, tol_norm):
    # Test hard-coded expected value without normalization
    assert edit.fast_birnbaum_simil(seq_x, seq_y) == expected

    # Test hard-coded expected value with normalization
    assert edit.fast_birnbaum_simil(seq_x, seq_y, normal=True) == pytest.approx(
        expected_norm, abs=tol_norm
    )


# TODO: add tests where fragile ends is different from levenshtein
# TODO: move fragile_end to similarity as it is not (necessarily) symmetric
# TODO: add a symmetric version with the mean?
@pytest.mark.parametrize(
    "seq_x,seq_y,expected,expected_norm,tol_norm",
    [
        ["kitten", "sitting"] + [3.0, 0.5, 0.0],
        [(1, 2, 3), [1, 2, 3]] + [0.0, 0.0, 0.0],
        [(1, 2, 3, 4, 5), (1, 2, 4, 3, 6, 7)] + [3.0, 0.6, 0.0],
        [(1, 2, 3), ["a", "b", "c", "d"]] + [4.0, 1.0, 0.0],
    ],
)
def test_fragile_ends_distance(seq_x, seq_y, expected, expected_norm, tol_norm):
    # Test hard-coded expected value without normalization
    assert edit.fragile_ends_simil(seq_x, seq_y) == pytest.approx(expected)

    # Test hard-coded expected value with normalization
    assert edit.fragile_ends_simil(seq_x, seq_y, normal=True) == pytest.approx(
        expected_norm, abs=tol_norm
    )

    # Test triangle-inequality
    seq_z = [element for element in seq_x] + [element for element in seq_y]
    dist_xz = edit.fragile_ends_simil(seq_x, seq_z)
    dist_xy = edit.fragile_ends_simil(seq_x, seq_y)
    dist_yz = edit.fragile_ends_simil(seq_y, seq_z)
    assert dist_xz <= (dist_xy + dist_yz)


# TODO: add tests where stemmatology is different from levenshtein
# TODO: move stemmatological to similarity as it is not (necessarily) symmetric
# TODO: add a symmetric version with the mean?
@pytest.mark.parametrize(
    "seq_x,seq_y,expected,expected_norm,tol_norm",
    [
        ["kitten", "sitting"] + [3.0, 0.428571, 1e-6],
        [(1, 2, 3), [1, 2, 3]] + [0.0, 0.0, 0.0],
        [(1, 2, 3, 4, 5), (1, 2, 4, 3, 6, 7)] + [3.0, 0.5, 0.0],
        [(1, 2, 3), ["a", "b", "c", "d"]] + [4.0, 1.0, 0.0],
    ],
)
def test_stemmatology_distance(seq_x, seq_y, expected, expected_norm, tol_norm):
    # Test hard-coded expected value without normalization
    assert edit.stemmatological_simil(seq_x, seq_y) == expected

    # Test hard-coded expected value without normalization
    assert edit.stemmatological_simil(seq_x, seq_y, normal=True) == pytest.approx(
        expected_norm, abs=tol_norm
    )

    # Test triangle-inequality
    seq_z = [element for element in seq_x] + [element for element in seq_y]
    dist_xz = edit.stemmatological_simil(seq_x, seq_z)
    dist_xy = edit.stemmatological_simil(seq_x, seq_y)
    dist_yz = edit.stemmatological_simil(seq_y, seq_z)
    assert dist_xz <= (dist_xy + dist_yz)


# TODO: add tests where stemmatology (norm/2030) is different from levenshtein (norm)
# TODO: move norm_stemmatological to similarity as it is not (necessarily) symmetric
# TODO: add a symmetric version with the mean?
# TODO: add non-normal test (if we keep this)
@pytest.mark.parametrize(
    "seq_x,seq_y,expected,tol",
    [
        ["kitten", "sitting"] + [0.428571, 1e-6],
        [(1, 2, 3), [1, 2, 3]] + [0.0, 0.0],
        [(1, 2, 3, 4, 5), (1, 2, 4, 3, 6, 7)] + [0.5, 0.0],
        [(1, 2, 3), ["a", "b", "c", "d"]] + [1.0, 0.0],
    ],
)
def test_stemmatology_2030_distance(seq_x, seq_y, expected, tol):
    # Test hard-coded expected value
    assert edit.stemmatological_simil(
        seq_x, seq_y, 20.0, 30.0, normal=True
    ) == pytest.approx(expected, abs=tol)

    # Test triangle-inequality
    seq_z = [element for element in seq_x] + [element for element in seq_y]
    dist_xz = edit.stemmatological_simil(seq_x, seq_z, 20.0, 30.0, normal=True)
    dist_xy = edit.stemmatological_simil(seq_x, seq_y, 20.0, 30.0, normal=True)
    dist_yz = edit.stemmatological_simil(seq_y, seq_z, 20.0, 30.0, normal=True)
    assert dist_xz <= (dist_xy + dist_yz)


# TESTS FOR MEASURES OF DISTANCE
# ------------------------------
#
# Measures of distance offer all the properties of distance measures (i.e.,


@pytest.mark.parametrize(
    "seq_x,seq_y,expected,expected_norm,tol_norm",
    [
        ["kitten", "sitting"] + [3.0, 0.428571, 1e-6],
        [(1, 2, 3), [1, 2, 3]] + [0.0, 0.0, 0.0],
        [(1, 2, 3, 4, 5), (1, 2, 4, 3, 6, 7)] + [3.0, 0.5, 0.0],
        [(1, 2, 3), ["a", "b", "c", "d"]] + [4.0, 1.0, 0.0],
    ],
)
def test_levenshtein_distance(seq_x, seq_y, expected, expected_norm, tol_norm):
    # Test hard-coded expected value without normalization
    assert edit.levenshtein_dist(seq_x, seq_y) == expected

    # Test hard-coded expected value with normalization
    assert edit.levenshtein_dist(seq_x, seq_y, normal=True) == pytest.approx(
        expected_norm, abs=tol_norm
    )

    # Test symmetry
    assert edit.levenshtein_dist(seq_x, seq_y) == edit.levenshtein_dist(seq_y, seq_x)

    # Test triangle-inequality
    seq_z = [element for element in seq_x] + [element for element in seq_y]
    dist_xz = edit.levenshtein_dist(seq_x, seq_z)
    dist_xy = edit.levenshtein_dist(seq_x, seq_y)
    dist_yz = edit.levenshtein_dist(seq_y, seq_z)
    assert dist_xz <= (dist_xy + dist_yz)


# TODO: add tests where levdamerau is different from levenshtein
@pytest.mark.parametrize(
    "seq_x,seq_y,expected,expected_norm,tol_norm",
    [
        ["kitten", "sitting"] + [3.0, 0.428571, 1e-6],
        [(1, 2, 3), [1, 2, 3]] + [0.0, 0.0, 0.0],
        [(1, 2, 3, 4, 5), (1, 2, 4, 3, 6, 7)] + [3.0, 0.5, 0.0],
        [(1, 2, 3), ["a", "b", "c", "d"]] + [4.0, 1.0, 0.0],
    ],
)
def test_levdamerau_distance(seq_x, seq_y, expected, expected_norm, tol_norm):
    # Test hard-coded expected value without normalization
    assert edit.levdamerau_dist(seq_x, seq_y) == expected

    # Test hard-coded expected value with normalization
    assert edit.levdamerau_dist(seq_x, seq_y, normal=True) == pytest.approx(
        expected_norm, abs=tol_norm
    )

    # Test symmetry
    assert edit.levdamerau_dist(seq_x, seq_y) == edit.levdamerau_dist(seq_y, seq_x)

    # Test triangle-inequality
    seq_z = [element for element in seq_x] + [element for element in seq_y]
    dist_xz = edit.levdamerau_dist(seq_x, seq_z)
    dist_xy = edit.levdamerau_dist(seq_x, seq_y)
    dist_yz = edit.levdamerau_dist(seq_y, seq_z)
    assert dist_xz <= (dist_xy + dist_yz)


# TODO: add tests where bulk delete is different from levenshtein
@pytest.mark.parametrize(
    "seq_x,seq_y,expected,expected_norm,tol_norm",
    [
        ["kitten", "sitting"] + [3.0, 0.428571, 1e-6],
        [(1, 2, 3), [1, 2, 3]] + [0.0, 0.0, 0.0],
        [(1, 2, 3, 4, 5), (1, 2, 4, 3, 6, 7)] + [3.0, 0.5, 0.0],
        [(1, 2, 3), ["a", "b", "c", "d"]] + [4.0, 1.0, 0.0],
    ],
)
def test_bulk_delete_distance(seq_x, seq_y, expected, expected_norm, tol_norm):
    # Test hard-coded expected value without normalization
    assert edit.bulk_delete_dist(seq_x, seq_y) == expected

    # Test hard-coded expected value without normalization
    assert edit.bulk_delete_dist(seq_x, seq_y, normal=True) == pytest.approx(
        expected_norm, abs=tol_norm
    )

    # Test symmetry
    assert edit.bulk_delete_dist(seq_x, seq_y) == edit.bulk_delete_dist(seq_y, seq_x)

    # Test triangle-inequality
    seq_z = [element for element in seq_x] + [element for element in seq_y]
    dist_xz = edit.bulk_delete_dist(seq_x, seq_z)
    dist_xy = edit.bulk_delete_dist(seq_x, seq_y)
    dist_yz = edit.bulk_delete_dist(seq_y, seq_z)
    assert dist_xz <= (dist_xy + dist_yz)


@pytest.mark.parametrize(
    "seq_x,seq_y,expected,tol",
    [
        ["kitten", "sitting"] + [0.253968, 1e-6],
        [(1, 2, 3), [1, 2, 3]] + [0.0, 0.0],
        [(1, 2, 3, 4, 5), (1, 2, 4, 3, 6, 7)] + [0.261111, 1e-6],
        [(1, 2, 3), ["a", "b", "c", "d"]] + [1.0, 0.0],
    ],
)
def test_jaro_distance(seq_x, seq_y, expected, tol):
    # Test hard-coded expected value
    assert edit.jaro_dist(seq_x, seq_y) == pytest.approx(expected, abs=tol)

    # Test symmetry
    assert edit.jaro_dist(seq_x, seq_y) == edit.jaro_dist(seq_y, seq_x)

    # Test triangle-inequality
    seq_z = [element for element in seq_x] + [element for element in seq_y]
    dist_xz = edit.jaro_dist(seq_x, seq_z)
    dist_xy = edit.jaro_dist(seq_x, seq_y)
    dist_yz = edit.jaro_dist(seq_y, seq_z)
    assert dist_xz <= (dist_xy + dist_yz)


@pytest.mark.parametrize(
    "seq_x,seq_y,expected,tol",
    [
        ["kitten", "sitting"] + [0.253968, 1e-6],
        [(1, 2, 3), [1, 2, 3]] + [0.0, 0.0],
        [(1, 2, 3, 4, 5), (1, 2, 4, 3, 6, 7)] + [0.208888, 1e-6],
        [(1, 2, 3), ["a", "b", "c", "d"]] + [1.0, 0.0],
    ],
)
def test_jarowinkler_distance(seq_x, seq_y, expected, tol):
    # Test hard-coded expected value
    assert edit.jaro_winkler_dist(seq_x, seq_y) == pytest.approx(expected, abs=tol)

    # Test symmetry
    assert edit.jaro_winkler_dist(seq_x, seq_y) == edit.jaro_winkler_dist(seq_y, seq_x)

    # Test triangle-inequality
    seq_z = [element for element in seq_x] + [element for element in seq_y]
    dist_xz = edit.jaro_winkler_dist(seq_x, seq_z)
    dist_xy = edit.jaro_winkler_dist(seq_x, seq_y)
    dist_yz = edit.jaro_winkler_dist(seq_y, seq_z)
    assert dist_xz <= (dist_xy + dist_yz)


@pytest.mark.parametrize(
    "seq_x,seq_y,expected,tol",
    [
        ["kitten", "sitting"] + [0.538461, 1e-6],
        [(1, 2, 3), [1, 2, 3]] + [0.0, 0.0],
        [(1, 2, 3, 4, 5), (1, 2, 4, 3, 6, 7)] + [0.554638, 1e-6],
        [(1, 2, 3), ["a", "b", "c", "d"]] + [1.0, 0.0],
    ],
)
def test_mmcwpa_distance(seq_x, seq_y, expected, tol):
    # Test hard-coded expected value
    assert edit.mmcwpa_dist(seq_x, seq_y) == pytest.approx(expected, abs=tol)

    # Test symmetry
    assert edit.mmcwpa_dist(seq_x, seq_y) == edit.mmcwpa_dist(seq_y, seq_x)

    # Test triangle-inequality
    seq_z = [element for element in seq_x] + [element for element in seq_y]
    dist_xz = edit.mmcwpa_dist(seq_x, seq_z)
    dist_xy = edit.mmcwpa_dist(seq_x, seq_y)
    dist_yz = edit.mmcwpa_dist(seq_y, seq_z)
    assert dist_xz <= (dist_xy + dist_yz)


@pytest.mark.parametrize(
    "seq_x,seq_y,expected,tol",
    [
        ["kitten", "sitting"] + [0.565217, 1e-6],
        [(1, 2, 3), [1, 2, 3]] + [0.0, 0.0],
        [(1, 2, 3, 4, 5), (1, 2, 4, 3, 6, 7)] + [0.666666, 1e-6],
        [(1, 2, 3), ["a", "b", "c", "d"]] + [1.0, 0.0],
    ],
)
def test_birnbaum_distance(seq_x, seq_y, expected, tol):
    # Test hard-coded expected value
    assert edit.birnbaum_dist(seq_x, seq_y) == pytest.approx(expected, abs=tol)

    # Test symmetry
    assert edit.birnbaum_dist(seq_x, seq_y) == edit.birnbaum_dist(seq_y, seq_x)

    # Test triangle-inequality
    seq_z = [element for element in seq_x] + [element for element in seq_y]
    dist_xz = edit.birnbaum_dist(seq_x, seq_z)
    dist_xy = edit.birnbaum_dist(seq_x, seq_y)
    dist_yz = edit.birnbaum_dist(seq_y, seq_z)
    assert dist_xz <= (dist_xy + dist_yz)


@pytest.mark.parametrize(
    "seq_x,seq_y,expected,tol",
    [
        ["kitten", "sitting"] + [0.666666, 1e-6],
        [(1, 2, 3), [1, 2, 3]] + [0.0, 0.0],
        [(1, 2, 3, 4, 5), (1, 2, 4, 3, 6, 7)] + [0.733333, 1e-6],
        [(1, 2, 3), ["a", "b", "c", "d"]] + [1.0, 0.0],
    ],
)
def test_fast_birnbaum_distance(seq_x, seq_y, expected, tol):
    # Test hard-coded expected value
    assert edit.fast_birnbaum_dist(seq_x, seq_y) == pytest.approx(expected, abs=tol)

    # Test symmetry
    assert edit.fast_birnbaum_dist(seq_x, seq_y) == edit.fast_birnbaum_dist(
        seq_y, seq_x
    )

    # Test triangle-inequality
    seq_z = [element for element in seq_x] + [element for element in seq_y]
    dist_xz = edit.fast_birnbaum_dist(seq_x, seq_z)
    dist_xy = edit.fast_birnbaum_dist(seq_x, seq_y)
    dist_yz = edit.fast_birnbaum_dist(seq_y, seq_z)
    assert dist_xz <= (dist_xy + dist_yz)
