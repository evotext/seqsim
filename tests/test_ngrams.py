"""
test_ngrams
===========

Tests for the `ngrams` module of the `seqsim` package.
"""

# Import Python standard libraries
import pickle

# Import the library being tested
from seqsim import ngrams
from seqsim.ngrams import PAD


def test_ngrams_iter_padded():
    assert list(ngrams.ngrams_iter("abc", 2)) == [
        (PAD, "a"),
        ("a", "b"),
        ("b", "c"),
        ("c", PAD),
    ]


def test_ngrams_iter_unpadded():
    assert list(ngrams.ngrams_iter("abcd", 3, pad_symbol=None)) == [
        ("a", "b", "c"),
        ("b", "c", "d"),
    ]


def test_ngrams_iter_custom_pad():
    assert list(ngrams.ngrams_iter([1, 2], 2, pad_symbol=0)) == [
        (0, 1),
        (1, 2),
        (2, 0),
    ]


def test_ngrams_pad_does_not_collide():
    # A sequence containing the old "$$$" pad string must be distinguishable
    # from padding
    grams = list(ngrams.ngrams_iter(["$$$", "a"], 2))
    assert grams == [(PAD, "$$$"), ("$$$", "a"), ("a", PAD)]
    assert PAD != "$$$"


def test_pad_is_singleton():
    assert pickle.loads(pickle.dumps(PAD)) is PAD
    assert repr(PAD) == "PAD"


def test_get_all_ngrams_by_order():
    assert list(ngrams.get_all_ngrams_by_order("ab", [1, 2], pad_symbol=None)) == [
        ("a",),
        ("b",),
        ("a", "b"),
    ]

    # Default collects all orders up to the sequence length
    all_grams = list(ngrams.get_all_ngrams_by_order("ab", pad_symbol=None))
    assert all_grams == [("a",), ("b",), ("a", "b")]
