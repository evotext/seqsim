"""
test_sequence
=============

Tests for the `sequence` module of the `seqsim` package.
"""

# Import Python standard libraries
import pytest

# Import the library being tested
from seqsim import sequence


@pytest.mark.parametrize(
    "seq_x,seq_y,expected,tol",
    [
        ["kitten", "sitting", 0.384615, 1e-6],
        [(1, 2, 3), [1, 2, 3], 0.0, 0.0],
        [(1, 2, 3, 4, 5), (1, 2, 4, 3, 6, 7), 0.454545, 1e-6],
        [(1, 2, 3), ["a", "b", "c", "d"], 1.0, 0.0],
    ],
)
def test_ratcliff_obershelp(seq_x, seq_y, expected, tol):
    # Test hard-coded expected value
    assert sequence.ratcliff_obershelp(seq_x, seq_y) == pytest.approx(expected, abs=tol)

    # Test symmetry
    assert sequence.ratcliff_obershelp(seq_x, seq_y) == sequence.ratcliff_obershelp(
        seq_y, seq_x
    )

    # Test triangle-inequality
    seq_z = [element for element in seq_x] + [element for element in seq_y]
    dist_xz = sequence.ratcliff_obershelp(seq_x, seq_z)
    dist_xy = sequence.ratcliff_obershelp(seq_x, seq_y)
    dist_yz = sequence.ratcliff_obershelp(seq_y, seq_z)
    assert dist_xz <= (dist_xy + dist_yz)
