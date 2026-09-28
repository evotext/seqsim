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
    "seq_x,seq_y,expected",
    [
        ("kitten", "sitting", 1 - 3 / 7),  # {i, n, t} / {e, g, i, k, n, s, t}
        ((1, 2, 3), [1, 2, 3], 0.0),
        ((1, 2, 3, 4, 5), (1, 2, 4, 3, 6, 7), 1 - 4 / 7),
        ((1, 2, 3), ["a", "b", "c", "d"], 1.0),
        ("ab", "ba", 0.0),  # order is ignored
    ],
)
def test_jaccard(seq_x, seq_y, expected):
    assert token.jaccard_dissim(seq_x, seq_y) == pytest.approx(expected)


@pytest.mark.parametrize(
    "seq_x,seq_y,expected",
    [
        # Multisets: kitten {k, i, t, t, e, n}, sitting {s, i, t, t, i, n, g};
        # intersection {i, t, t, n}
        ("kitten", "sitting", 1 - 8 / 13),
        ((1, 2, 3), [1, 2, 3], 0.0),
        ((1, 2, 3, 4, 5), (1, 2, 4, 3, 6, 7), 1 - 8 / 11),
        ((1, 2, 3), ["a", "b", "c", "d"], 1.0),
        ("aab", "abb", 1 - 4 / 6),
    ],
)
def test_sorensen(seq_x, seq_y, expected):
    assert token.sorensen_dissim(seq_x, seq_y) == pytest.approx(expected)


@pytest.mark.parametrize(
    "seq_x,seq_y,expected",
    [
        # Length 1: 2/5, length 2: 1/4, lengths 3 and 4: 0; weights 1..4
        ("abc", "bcde", 1 - (0.4 * 1 + 0.25 * 2) / 10),
        ((1, 2, 3), [1, 2, 3], 0.0),
        ((1, 2, 3), ["a", "b", "c", "d"], 1.0),
        # Length 1: {a, a} vs {a} gives 1/2, length 2: {aa} vs {} gives 0
        ("aa", "a", 1 - (0.5 * 1) / 3),
    ],
)
def test_subseq_jaccard(seq_x, seq_y, expected):
    assert token.subseq_jaccard_dissim(seq_x, seq_y) == pytest.approx(expected)


@pytest.mark.parametrize("seq", ["aaa", "abab", "abcabc", (1, 1, 2, 1, 1)])
def test_subseq_jaccard_identity(seq):
    # In 0.3.1, repeated sub-sequences gave a non-zero score for identical
    # sequences
    assert token.subseq_jaccard_dissim(seq, seq) == 0.0
