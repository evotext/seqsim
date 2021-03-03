"""
test_general
============

Tests for the common wrappers offered by the `seqsim` package. The
other individual tests (e.g., `test_edit-py`) are designed for a
more detailed testing of the individual methods, mostly working with
short sequences, almost always in pairwise, and with a focus on the
results. The tests in this module are designed more for coverage,
profiling, and bug-tracking, also as a very simplified fuzzer, while also
testing multiple sequences.
"""

# Import Python standard libraries
import pytest

# Import the library being tested
import seqsim


@pytest.mark.parametrize(
    "seq_x,seq_y",
    [
        ["test", "test"],
        ["test", "tset"],
        ["test", "testest"],
        ["test", "testtest"],
        ["aaa", "bbb"],
        ["cat", "hat"],
        ["Niall", "Neil"],
        ["aluminum", "Catalan"],
        ["ATCG", "TAGC"],
        ["GATTACA", "GCATGCU"],
        ["AGACTAGTTAC", "TGACGSTGC"],
        ["AGACTAGTTAC", "CGAGACGT"],
    ],
)
def test_pairwise_distance(seq_x, seq_y):
    """
    Test all methods pairwise.
    """

    for method in seqsim.METHODS:
        assert seqsim.distance([seq_x, seq_y], method=method) >= 0
        assert seqsim.distance([seq_x, seq_y], method=method, normal=True) >= 0


@pytest.mark.parametrize(
    "seq_x,seq_y,seq_z",
    [
        ["test", "test", "test"],
        ["test", "tset", "test"],
        ["test", "testest", "ttest"],
        ["test", "testtest", "testtesttest"],
        ["aaa", "bbb", "ccc"],
        ["cat", "hat", "hac"],
        ["Niall", "Neil", "neil"],
        ["aluminum", "Catalan", "Italian"],
        ["ATCG", "TAGC", "CGAT"],
        ["GATTACA", "GCATGCU", "GATAC"],
        ["AGACTAGTTAC", "TGACGSTGC", "TGACGTGCAATTA"],
        ["AGACTAGTTAC", "CGAGACGT", "CGAGACTTTTTTTTTTT"],
    ],
)
def test_multiwise_distance(seq_x, seq_y, seq_z):
    """
    Test all methods multiwise.
    """

    for method in seqsim.METHODS:
        assert seqsim.distance([seq_x, seq_y, seq_z], method=method) >= 0
        assert seqsim.distance([seq_x, seq_y, seq_z], method=method, normal=True) >= 0
