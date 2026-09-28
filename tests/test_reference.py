"""
test_reference
==============

Tests comparing the optimized implementations of some methods against
straightforward reference implementations, which follow the published
descriptions of the algorithms literally.
"""

# Import Python standard libraries
from collections import Counter

import pytest
from hypothesis import given, settings, strategies as st

# Import the library being tested
from seqsim import edit, token

elements = st.sampled_from(["a", "b", "c", 1, (2, 3)])
sequences = st.lists(elements, max_size=14)


def _reference_mmcwpa_ssnc(seq_x, seq_y):
    """
    MMCWPA as in the original description, trying all window lengths in
    decreasing order and searching sub-sequences element by element.
    """

    def find(hay, needle):
        for i in range(len(hay) - len(needle) + 1):
            if tuple(hay[i : i + len(needle)]) == tuple(needle):
                return i
        return None

    f_x, f_y, ssnc = [list(seq_x)], [list(seq_y)], 0
    while f_x and f_y:
        found = None
        for idx_x, sf_x in enumerate(f_x):
            for length in range(len(sf_x), 0, -1):
                for i in range(len(sf_x) - length + 1):
                    for idx_y, sf_y in enumerate(f_y):
                        j = find(sf_y, sf_x[i : i + length])
                        if j is not None:
                            found = (idx_x, i, idx_y, j, length)
                            break
                    if found:
                        break
                if found:
                    break
            if found:
                break
        if not found:
            break

        idx_x, i, idx_y, j, length = found
        sf_x, sf_y = f_x[idx_x], f_y[idx_y]
        f_x[idx_x : idx_x + 1] = [s for s in (sf_x[:i], sf_x[i + length :]) if s]
        f_y[idx_y : idx_y + 1] = [s for s in (sf_y[:j], sf_y[j + length :]) if s]
        ssnc += (2 * length) ** 2

    return ssnc


def _reference_subseq_jaccard(seq_x, seq_y):
    """
    Subseq-Jaccard computed directly on tuples.
    """

    max_length = max(len(seq_x), len(seq_y))
    total = 0.0
    for length in range(1, max_length + 1):
        grams_x = Counter(
            tuple(seq_x[i : i + length]) for i in range(len(seq_x) - length + 1)
        )
        grams_y = Counter(
            tuple(seq_y[i : i + length]) for i in range(len(seq_y) - length + 1)
        )
        inter = sum((grams_x & grams_y).values())
        union = sum((grams_x | grams_y).values())
        total += length * inter / union

    return 1.0 - total / (max_length * (max_length + 1) / 2)


@given(seq_x=sequences.filter(len), seq_y=sequences.filter(len))
@settings(max_examples=300, deadline=None)
def test_mmcwpa_matches_reference(seq_x, seq_y):
    ssnc = max(
        _reference_mmcwpa_ssnc(seq_x, seq_y), _reference_mmcwpa_ssnc(seq_y, seq_x)
    )
    expected = 1.0 - (ssnc**0.5) / (len(seq_x) + len(seq_y))
    assert edit.mmcwpa_dissim(seq_x, seq_y) == pytest.approx(expected, abs=1e-12)


@given(seq_x=sequences.filter(len), seq_y=sequences.filter(len))
@settings(max_examples=300, deadline=None)
def test_subseq_jaccard_matches_reference(seq_x, seq_y):
    assert token.subseq_jaccard_dissim(seq_x, seq_y) == pytest.approx(
        _reference_subseq_jaccard(seq_x, seq_y), abs=1e-12
    )
