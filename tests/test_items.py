"""
test_items
==========

Tests for the views of a sequence as items (`seqsim._items`).
"""

# Import Python standard libraries
import pickle
from collections import Counter

import pytest

# Import the library being tested
from seqsim import _items, ngrams
from seqsim._items import BOUNDARY


def test_occurrences():
    assert _items.occurrences("abab") == [("a", 0), ("b", 0), ("a", 1), ("b", 1)]
    assert _items.occurrences([]) == []


def test_first_occurrences():
    assert _items.first_occurrences("abacb") == ["a", "b", "c"]


@pytest.mark.parametrize(
    "repeats,expected",
    [
        (
            "occurrence",
            ([("a", 0), ("b", 0), ("b", 1)], [("b", 0), ("a", 0), ("b", 1)]),
        ),
        ("first", ([("a", 0), ("b", 0)], [("b", 0), ("a", 0)])),
    ],
)
def test_shared(repeats, expected):
    assert _items.shared("abab", "bab", repeats=repeats) == expected


def test_shared_invalid_repeats():
    with pytest.raises(ValueError):
        _items.shared("ab", "ab", repeats="last")


def test_adjacencies():
    assert _items.adjacencies("aba", boundaries=False) == Counter(
        {("a", "b"): 1, ("b", "a"): 1}
    )
    assert _items.adjacencies("ab", boundaries=True) == Counter(
        {(BOUNDARY, "a"): 1, ("a", "b"): 1, ("b", BOUNDARY): 1}
    )
    # A sequence of `n` elements has `n + 1` adjacencies with the boundaries
    assert sum(_items.adjacencies("", boundaries=True).values()) == 1
    assert sum(_items.adjacencies("abcab", boundaries=True).values()) == 6


def test_ngrams():
    assert _items.ngrams("abab", 2) == Counter({("a", "b"): 2, ("b", "a"): 1})
    assert _items.ngrams("a", 2) == Counter()
    assert _items.ngrams("a", 2, pad=True) == Counter(
        {(BOUNDARY, "a"): 1, ("a", BOUNDARY): 1}
    )
    assert _items.ngrams("abc", 3, pad=True) == Counter(ngrams.ngrams_iter("abc", 3))


def test_boundary_is_a_single_sentinel():
    assert BOUNDARY is ngrams.PAD
    assert pickle.loads(pickle.dumps(BOUNDARY)) is BOUNDARY
    assert BOUNDARY != "PAD"
