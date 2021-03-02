"""
test_token
==========

Tests for the `token` module of the `seqsim` package.
"""

# Import Python standard libraries
import pytest

# Import the library being tested
from seqsim import token


@pytest.mark.parametrize(
    "seq_x,seq_y,expected,tol",
    [
        ["kitten", "sitting"] + [0.7, 0.0],
        [(1, 2, 3), [1, 2, 3]] + [0.0, 0.0],
        [(1, 2, 3, 4, 5), (1, 2, 4, 3, 6, 7)] + [0.428571, 1e-6],
        [(1, 2, 3), ["a", "b", "c", "d"]] + [1.0, 0.0],
    ],
)
def test_jaccard_distance(seq_x, seq_y, expected, tol):
    # Test hard-coded expected value
    assert token.jaccard_dist(seq_x, seq_y) == pytest.approx(expected, abs=tol)

    # Test symmetry
    assert token.jaccard_dist(seq_x, seq_y) == token.jaccard_dist(seq_y, seq_x)

    # Test triangle-inequality
    seq_z = [element for element in seq_x] + [element for element in seq_y]
    dist_xz = token.jaccard_dist(seq_x, seq_z)
    dist_xy = token.jaccard_dist(seq_x, seq_y)
    dist_yz = token.jaccard_dist(seq_y, seq_z)
    assert dist_xz <= (dist_xy + dist_yz)


@pytest.mark.parametrize(
    "seq_x,seq_y,expected,tol",
    [
        ["kitten", "sitting"] + [0.751556, 1e-6],
        [(1, 2, 3), [1, 2, 3]] + [0.0, 0.0],
        [(1, 2, 3, 4, 5), (1, 2, 4, 3, 6, 7)] + [0.787094, 1e-6],
        [(1, 2, 3), ["a", "b", "c", "d"]] + [1.0, 0.0],
    ],
)
def test_subseq_jaccard_distance(seq_x, seq_y, expected, tol):
    # Test hard-coded expected value
    assert token.subseq_jaccard_dist(seq_x, seq_y) == pytest.approx(expected, abs=tol)

    # Test symmetry
    assert token.subseq_jaccard_dist(seq_x, seq_y) == token.subseq_jaccard_dist(
        seq_y, seq_x
    )

    # Test triangle-inequality
    seq_z = [element for element in seq_x] + [element for element in seq_y]
    dist_xz = token.subseq_jaccard_dist(seq_x, seq_z)
    dist_xy = token.subseq_jaccard_dist(seq_x, seq_y)
    dist_yz = token.subseq_jaccard_dist(seq_y, seq_z)
    assert dist_xz <= (dist_xy + dist_yz)
