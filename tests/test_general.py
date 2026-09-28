"""
test_general
============

Tests for the common wrappers offered by the `seqsim` package. The
other individual tests (e.g., `test_edit.py`) are designed for a
more detailed testing of the individual methods, mostly working with
short sequences, almost always in pairwise, and with a focus on the
results. The tests in this module are designed more for coverage and
the `distance()` wrapper, also testing multiple sequences.
"""

# Import Python standard libraries
import itertools

import pytest

# Import the library being tested
import seqsim

WORDS = [
    "test",
    "tset",
    "testest",
    "testtesttest",
    "aaa",
    "bbb",
    "cat",
    "hat",
    "Niall",
    "Neil",
    "aluminum",
    "Catalan",
    "ATCG",
    "TAGC",
    "GATTACA",
    "GCATGCU",
    "AGACTAGTTAC",
]


@pytest.mark.parametrize("method", sorted(seqsim.METHODS))
def test_pairwise_distance(method):
    for seq_x, seq_y in itertools.combinations(WORDS, 2):
        dist = seqsim.distance([seq_x, seq_y], method=method)
        dist_norm = seqsim.distance([seq_x, seq_y], method=method, normal=True)
        assert isinstance(dist, float)
        assert dist >= 0.0
        assert 0.0 <= dist_norm <= 1.0


@pytest.mark.parametrize("method", sorted(seqsim.METHODS))
def test_multiwise_distance(method):
    for seqs in itertools.combinations(WORDS[:8], 3):
        pairwise = [
            seqsim.distance([seq_x, seq_y], method=method)
            for seq_x, seq_y in itertools.combinations(seqs, 2)
        ]
        assert seqsim.distance(seqs, method=method) == pytest.approx(
            sum(pairwise) / len(pairwise)
        )


def test_distance_many_sequences():
    # The mean is taken over all pairs, also for more than three sequences
    assert seqsim.distance(["abc", "abc", "abc", "xyz"]) == pytest.approx(1.5)
    assert seqsim.distance(list("abcde"), normal=True) == pytest.approx(1.0)


def test_distance_iterables():
    # Generators, of sequences and inside sequences, are accepted
    seqs = (seq for seq in ["kitten", "sitting"])
    assert seqsim.distance(seqs) == 3.0
    assert seqsim.distance([iter("kitten"), iter("sitting")]) == 3.0


def test_distance_kwargs():
    seqs = ["abcdeXXXXXfghij", "abcdefghij"]
    assert seqsim.distance(seqs, "bulk_delete") == 1.0
    assert seqsim.distance(seqs, "bulk_delete", max_del_len=1) == 5.0
    assert (
        seqsim.distance(
            seqs, "stemmatological", frag_start=0.0, frag_end=0.0, max_del_len=2
        )
        == 3.0
    )


def test_distance_errors():
    with pytest.raises(TypeError):
        seqsim.distance("ab")
    with pytest.raises(ValueError):
        seqsim.distance(["ab"])
    with pytest.raises(ValueError):
        seqsim.distance(["ab", "cd"], method="unknown")
    with pytest.raises(TypeError):
        seqsim.distance(["ab", "cd"], method="levenshtein", max_del_len=2)


@pytest.mark.parametrize("method", sorted(seqsim.METHODS))
def test_keyword_only_parameters(method):
    # Parameters other than the sequences must be passed by name
    with pytest.raises(TypeError):
        seqsim.METHODS[method]("abc", "abd", True)
