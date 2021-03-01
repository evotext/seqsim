"""
test_similarity
===============

Tests for the `similarity` module of the `seqsim` package.
"""

# Import Python standard libraries
import pytest

# Import the library being tested
from seqsim import similarity

TEST1 = ["kitten", "sitting"]
TEST2 = [(1, 2, 3), [1, 2, 3]]
TEST3 = [(1, 2, 3, 4, 5), (1, 2, 4, 3, 6, 7)]
TEST4 = [(1, 2, 3), ["a", "b", "c", "d"]]


@pytest.mark.parametrize(
    "seq_x,seq_y,expected",
    [
        TEST1 + [7],
        TEST2 + [6],
        TEST3 + [4],
        TEST4 + [0],
    ],
)
def test_fast_birnbaum_simil(seq_x, seq_y, expected):
    # Test hard-coded expected value
    assert similarity.fast_birnbaum_simil(seq_x, seq_y) == expected


@pytest.mark.parametrize(
    "seq_x,seq_y,expected",
    [
        TEST1 + [10],
        TEST2 + [6],
        TEST3 + [5],
        TEST4 + [0],
    ],
)
def test_birnbaum_simil(seq_x, seq_y, expected):
    # Test hard-coded expected value
    assert similarity.birnbaum_simil(seq_x, seq_y) == expected
