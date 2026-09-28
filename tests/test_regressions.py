"""
test_regressions
================

Regression tests for bugs found in the adversarial review of version 0.3.1.
Each test fails on 0.3.1 (with the functions under their old names).
"""

# Import Python standard libraries
import pytest

# Import the library being tested
import seqsim
from seqsim import edit, sequence, token


def test_mmcwpa_searches_all_subfields():
    # 0.3.1 stopped after the first subfield of Fx without a match
    assert edit.mmcwpa_dissim("QabRcd", "abcd") == pytest.approx(0.434315, abs=1e-6)
    assert edit.mmcwpa_dissim("abcd", "QabRcd") == pytest.approx(0.434315, abs=1e-6)
    assert edit.mmcwpa_dissim("kitten", "sitting") == pytest.approx(0.513496, abs=1e-6)


def test_distance_mean_over_pairs():
    # 0.3.1 divided the sum of pairwise distances by the number of sequences
    assert seqsim.distance(["abc", "abc", "abc", "xyz"]) == pytest.approx(1.5)
    assert seqsim.distance(list("abcde"), normal=True) == pytest.approx(1.0)


def test_stemmatological_full_block_deletion():
    # A block of exactly `max_del_len` elements is a single deletion
    assert edit.stemmatological_dissim(
        "abcdeXXXXXfghij", "abcdefghij", frag_start=0.0, frag_end=0.0
    ) == pytest.approx(1.0)

    # With `max_del_len=1`, deletions behave as in Levenshtein
    assert edit.stemmatological_dissim(
        "abcXdef", "abcdef", frag_start=0.0, frag_end=0.0, max_del_len=1
    ) == pytest.approx(1.0)


@pytest.mark.parametrize("max_del_len", [0, -1])
def test_invalid_max_del_len(max_del_len):
    with pytest.raises(ValueError):
        edit.bulk_delete_dist("abc", "abd", max_del_len=max_del_len)


def test_birnbaum_simil_identical_normalized():
    assert edit.birnbaum_simil("abc", "abc", normal=True) == 1.0
    assert edit.birnbaum_simil((1, 2, 3), [1, 2, 3], normal=True) == 1.0


def test_ratcliff_obershelp_no_str_collisions():
    assert sequence.ratcliff_obershelp_dissim([1, 23], [12, 3]) == 1.0


def test_subseq_jaccard_identity():
    assert token.subseq_jaccard_dissim("aaa", "aaa") == 0.0


@pytest.mark.parametrize("method", sorted(seqsim.METHODS))
def test_empty_sequences_do_not_crash(method):
    # Many methods raised ZeroDivisionError on empty sequences
    assert seqsim.distance(["", ""], method, normal=True) == 0.0
    assert seqsim.distance(["", ""], method) == 0.0
    assert seqsim.distance(["a", ""], method, normal=True) > 0.0


@pytest.mark.parametrize(
    "func,seq_x,seq_y",
    [
        # Asymmetric cases reported in the review
        (edit.bulk_delete_dist, "abcdefgh", "a"),
        (edit.birnbaum_dissim, "abcbc", "bbbba"),
        (edit.fragile_ends_dissim, "abcdefghij", "bcdefghij"),
        (edit.mmcwpa_dissim, "caab", "ba"),
    ],
)
def test_symmetry(func, seq_x, seq_y):
    assert func(seq_x, seq_y) == func(seq_y, seq_x)


def test_damerau_triangle_inequality():
    # 0.3.1 computed OSA, where d("ca", "abc") = 3 > d("ca", "ac") + d("ac", "abc")
    assert edit.damerau_dist("ca", "abc") <= edit.damerau_dist(
        "ca", "ac"
    ) + edit.damerau_dist("ac", "abc")
