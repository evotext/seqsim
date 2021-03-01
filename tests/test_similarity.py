"""
test_similarity
===============

Tests for the `similarity` module of the `seqsim` package.
"""

# Import Python standard libraries
import pytest

# Import the library being tested
from seqsim import similarity

test1_a = "kitten"
test1_b = "sitting"
test2_a = (1, 2, 3)
test2_b = [1, 2, 3]
test3_a = (1, 2, 3, 4, 5)
test3_b = (1, 2, 4, 3, 6, 7)
test4_a = (1, 2, 3)
test4_b = ["a", "b", "c", "d"]


def test_birnbaum_similarity():
    assert similarity.fast_birnbaum_simil(test1_a, test1_b) == 7
    assert similarity.fast_birnbaum_simil(test2_a, test2_b) == 6
    assert similarity.fast_birnbaum_simil(test3_a, test3_b) == 4
    assert similarity.fast_birnbaum_simil(test4_a, test4_b) == 0
