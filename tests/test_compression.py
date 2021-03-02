"""
test_compression
================

Tests for the `compression` module of the `seqsim` package.
"""

# Import Python standard libraries
import pytest

# Import the library being tested
from seqsim import compression


@pytest.mark.parametrize(
    "seq_x,seq_y,expected,tol",
    [
        ["kitten", "sitting", 1.25, 0.0],
        [(1, 2, 3), [1, 2, 3], 1.25, 0.0],
        [(1, 2, 3, 4, 5), (1, 2, 4, 3, 6, 7), 1.0, 0.0],
        [(1, 2, 3), ["a", "b", "c", "d"], 1.555555, 1e-6],
    ],
)
def test_arith_ncd(seq_x, seq_y, expected, tol):
    # Test hard-coded expected value
    assert compression.arith_ncd(seq_x, seq_y) == pytest.approx(expected, abs=tol)

    # Test symmetry
    assert compression.arith_ncd(seq_x, seq_y) == compression.arith_ncd(seq_y, seq_x)


@pytest.mark.parametrize(
    "seq_x,seq_y,expected,tol",
    [
        ["kitten", "sitting", 0.101340, 1e-6],
        [(1, 2, 3), [1, 2, 3], 0.0, 0.0],
        [(1, 2, 3, 4, 5), (1, 2, 4, 3, 6, 7), 0.114430, 1e-6],
        [(1, 2, 3), ["a", "b", "c", "d"], 0.407464, 1e-6],
    ],
)
def test_entropy_ncd(seq_x, seq_y, expected, tol):
    # Test hard-coded expected value
    assert compression.entropy_ncd(seq_x, seq_y) == pytest.approx(expected, abs=tol)

    # Test symmetry
    assert compression.entropy_ncd(seq_x, seq_y) == compression.entropy_ncd(
        seq_y, seq_x
    )
