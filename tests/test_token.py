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


@pytest.mark.parametrize(
    "seq_x,seq_y,kwargs,expected",
    [
        ("abc", "bca", {}, 6.0),
        ("abc", "bca", {"pad": False}, 2.0),
        ("abaca", "acaba", {}, 0.0),  # same profile, different sequences
        ("abcde", "abcde", {"q": 3}, 0.0),
        ("a", "", {}, 3.0),  # (PAD, a) and (a, PAD) vs (PAD, PAD)
        ("a", "", {"pad": False}, 0.0),  # no 2-grams without padding
        ("abc", "xyz", {"q": 1, "pad": False}, 6.0),
    ],
)
def test_qgram(seq_x, seq_y, kwargs, expected):
    assert token.qgram_dissim(seq_x, seq_y, **kwargs) == expected
    assert token.qgram_dissim(seq_y, seq_x, **kwargs) == expected


def test_qgram_normal_and_errors():
    assert token.qgram_dissim("abc", "xyz", normal=True) == 1.0
    assert token.qgram_dissim("", "", normal=True) == 0.0
    with pytest.raises(ValueError):
        token.qgram_dissim("abc", "abd", q=0)


@pytest.mark.parametrize(
    "seq_x,seq_y,alpha,beta,expected",
    [
        ("abc", "abcdef", 0.5, 0.5, 2 * 3 / 9),  # Sørensen–Dice
        ("abc", "abcdef", 1.0, 1.0, 3 / 6),  # Jaccard
        ("abc", "abcdef", 1.0, 0.0, 1.0),  # all of x is in y
        ("abcdef", "abc", 1.0, 0.0, 0.5),
        ("", "", 0.5, 0.5, 1.0),
        ("abc", "", 0.5, 0.5, 0.0),
    ],
)
def test_tversky(seq_x, seq_y, alpha, beta, expected):
    assert token.tversky_simil(seq_x, seq_y, alpha=alpha, beta=beta) == pytest.approx(
        expected
    )


def test_tversky_matches_sorensen():
    for seq_x, seq_y in [("kitten", "sitting"), ("aab", "abb"), ("abc", "xyz")]:
        assert token.tversky_simil(seq_x, seq_y) == pytest.approx(
            1 - token.sorensen_dissim(seq_x, seq_y)
        )
    with pytest.raises(ValueError):
        token.tversky_simil("abc", "abd", alpha=-1.0)


@pytest.mark.parametrize(
    "seq_x,seq_y,size,expected",
    [
        ("bcd", "abcdef", 2, 1.0),
        ("abcdef", "bcd", 2, 0.4),
        ("abc", "cba", 1, 1.0),
        ("abc", "cba", 2, 0.0),
        ("", "abc", 1, 1.0),
        ("abc", "", 1, 0.0),
    ],
)
def test_containment(seq_x, seq_y, size, expected):
    assert token.containment(seq_x, seq_y, size=size) == pytest.approx(expected)


def test_qgram_ukkonen_example():
    # Ukkonen (1992), p. 193
    assert token.qgram_dissim("01000", "001111", q=2, pad=False) == 5.0
