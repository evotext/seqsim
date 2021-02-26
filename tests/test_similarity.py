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
    assert similarity.birnbaum(test1_a, test1_b) == 8
    assert similarity.birnbaum(test2_a, test2_b) == 7
    assert similarity.birnbaum(test3_a, test3_b) == 5
    assert similarity.birnbaum(test4_a, test4_b) == 1
