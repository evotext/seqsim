"""
test_properties
===============

Property-based tests (with `hypothesis`) of the mathematical properties that
the name of each method promises:

  * all methods: non-negativity, `0.0` for identical sequences, normalized
    values in range [0..1], float return values, and defined values for
    empty sequences;
  * all `_dist` and `_dissim` methods: symmetry;
  * `_dist` methods: identity of indiscernibles and the triangle inequality;
  * `_simil` methods: symmetry and maximum (normalized) similarity only for
    identical sequences.
"""

# Import Python standard libraries
import math

import pytest
from hypothesis import given, settings, strategies as st

# Import the library being tested
from seqsim import alignment, compression, edit, order, sequence, token

DISTANCES = [
    edit.indel_dist,
    edit.lcs_dist,
    edit.levenshtein_gld_dist,
    edit.damerau_gld_dist,
    edit.indel_gld_dist,
    edit.levenshtein_ned_dist,
    order.ulam_dist,
    edit.levenshtein_dist,
    edit.damerau_dist,
    edit.bulk_delete_dist,
]

DISSIMILARITIES = [
    edit.block_move_dissim,
    edit.gst_dissim,
    order.kendall_tau_dissim,
    order.footrule_dissim,
    order.cayley_dissim,
    order.block_interchange_dissim,
    order.breakpoint_dissim,
    token.qgram_dissim,
    compression.lz76_dissim,
    alignment.nw_dissim,
    edit.osa_dissim,
    edit.fragile_ends_dissim,
    edit.stemmatological_dissim,
    edit.jaro_dissim,
    edit.jaro_winkler_dissim,
    edit.mmcwpa_dissim,
    edit.birnbaum_dissim,
    token.jaccard_dissim,
    token.sorensen_dissim,
    token.subseq_jaccard_dissim,
    sequence.ratcliff_obershelp_dissim,
    compression.entropy_ncd_dissim,
    compression.lzma_ncd_dissim,
]

# Dissimilarities for which identical sequences, and only identical
# sequences, score zero (LZMA NCD is positive for identical short sequences,
# and Jaccard, Sørensen and entropy NCD ignore order)
DISSIM_WITH_IDENTITY = [
    edit.block_move_dissim,
    order.kendall_tau_dissim,
    order.footrule_dissim,
    order.cayley_dissim,
    order.block_interchange_dissim,
    alignment.nw_dissim,
    edit.osa_dissim,
    edit.fragile_ends_dissim,
    edit.stemmatological_dissim,
    edit.jaro_dissim,
    edit.jaro_winkler_dissim,
    edit.mmcwpa_dissim,
    edit.birnbaum_dissim,
    token.subseq_jaccard_dissim,
    sequence.ratcliff_obershelp_dissim,
]

ALL = DISTANCES + DISSIMILARITIES

# Small alphabets make shared elements, and thus interesting cases, likely;
# elements of different types are mixed to exercise arbitrary hashables
elements = st.sampled_from(["a", "b", "c", 1, (2, 3), None])
seqs = st.lists(elements, max_size=8)
strings = st.text(alphabet="abc", max_size=8)
sequences = st.one_of(seqs, strings, seqs.map(tuple))

ids = lambda func: func.__name__  # noqa: E731


def close(value, target):
    return math.isclose(value, target, abs_tol=1e-9)


@pytest.mark.parametrize("func", ALL, ids=ids)
@given(seq_x=sequences, seq_y=sequences)
@settings(max_examples=150, deadline=None)
def test_range_and_type(func, seq_x, seq_y):
    value = func(seq_x, seq_y)
    value_norm = func(seq_x, seq_y, normal=True)

    assert isinstance(value, float)
    assert isinstance(value_norm, float)
    assert value >= 0.0
    assert 0.0 <= value_norm <= 1.0


@pytest.mark.parametrize("func", ALL, ids=ids)
@given(seq_x=sequences, seq_y=sequences)
@settings(max_examples=150, deadline=None)
def test_symmetry(func, seq_x, seq_y):
    assert close(func(seq_x, seq_y), func(seq_y, seq_x))
    assert close(func(seq_x, seq_y, normal=True), func(seq_y, seq_x, normal=True))


@pytest.mark.parametrize(
    "func",
    [f for f in ALL if f not in (compression.lzma_ncd_dissim, compression.lz76_dissim)],
    ids=ids,
)
@given(seq=sequences)
@settings(max_examples=100, deadline=None)
def test_identical_is_zero(func, seq):
    assert func(seq, list(seq)) == 0.0
    assert func(seq, list(seq), normal=True) == 0.0


@pytest.mark.parametrize("func", DISTANCES + DISSIM_WITH_IDENTITY, ids=ids)
@given(seq_x=sequences, seq_y=sequences)
@settings(max_examples=150, deadline=None)
def test_identity_of_indiscernibles(func, seq_x, seq_y):
    if list(seq_x) != list(seq_y):
        assert func(seq_x, seq_y) > 0.0
        assert func(seq_x, seq_y, normal=True) > 0.0


@pytest.mark.parametrize("func", DISTANCES, ids=ids)
@given(seq_x=sequences, seq_y=sequences, seq_z=sequences)
@settings(max_examples=300, deadline=None)
def test_triangle_inequality(func, seq_x, seq_y, seq_z):
    assert func(seq_x, seq_z) <= func(seq_x, seq_y) + func(seq_y, seq_z) + 1e-9


@pytest.mark.parametrize("max_del_len", [1, 2, 3])
@given(seq_x=strings, seq_y=strings, seq_z=strings)
@settings(max_examples=200, deadline=None)
def test_bulk_delete_triangle_inequality(max_del_len, seq_x, seq_y, seq_z):
    def dist(a, b):
        return edit.bulk_delete_dist(a, b, max_del_len=max_del_len)

    assert dist(seq_x, seq_z) <= dist(seq_x, seq_y) + dist(seq_y, seq_z)


@pytest.mark.parametrize("func", DISTANCES + DISSIMILARITIES, ids=ids)
@given(seq=st.one_of(seqs, strings).filter(len))
@settings(max_examples=50, deadline=None)
def test_empty_sequences(func, seq):
    empty = type(seq)()
    assert func(empty, empty) == 0.0
    assert func(empty, empty, normal=True) == 0.0
    assert func(seq, empty) > 0.0
    assert func(seq, empty, normal=True) > 0.0


@given(seq_x=sequences, seq_y=sequences)
@settings(max_examples=150, deadline=None)
def test_birnbaum_simil(seq_x, seq_y):
    simil = edit.birnbaum_simil(seq_x, seq_y)
    simil_norm = edit.birnbaum_simil(seq_x, seq_y, normal=True)

    assert isinstance(simil, float)
    assert simil == edit.birnbaum_simil(seq_y, seq_x)
    assert 0.0 <= simil_norm <= 1.0
    assert (simil_norm == 1.0) == (list(seq_x) == list(seq_y))


@given(seq_x=sequences, seq_y=sequences)
@settings(max_examples=150, deadline=None)
def test_sw_simil(seq_x, seq_y):
    simil = alignment.sw_simil(seq_x, seq_y)
    simil_norm = alignment.sw_simil(seq_x, seq_y, normal=True)

    assert isinstance(simil, float)
    assert simil == alignment.sw_simil(seq_y, seq_x)
    assert 0.0 <= simil_norm <= 1.0
    if list(seq_x) == list(seq_y):
        assert simil_norm == 1.0


@given(
    seq_x=st.lists(strings, max_size=4),
    seq_y=st.lists(strings, max_size=4),
)
@settings(max_examples=100, deadline=None)
def test_monge_elkan_simil(seq_x, seq_y):
    simil = alignment.monge_elkan_simil(seq_x, seq_y)

    assert 0.0 <= simil <= 1.0
    assert math.isclose(simil, alignment.monge_elkan_simil(seq_y, seq_x))
    if list(seq_x) == list(seq_y):
        assert simil == 1.0


def breakpoint_inner(seq_x, seq_y, **kwargs):
    return order.breakpoint_dissim(seq_x, seq_y, boundaries=False, **kwargs)


@given(seq_x=sequences, seq_y=sequences, seq_z=sequences)
@settings(max_examples=300, deadline=None)
def test_breakpoint_without_boundaries(seq_x, seq_y, seq_z):
    value = breakpoint_inner(seq_x, seq_y)
    value_norm = breakpoint_inner(seq_x, seq_y, normal=True)

    assert isinstance(value, float)
    assert value >= 0.0
    assert 0.0 <= value_norm <= 1.0
    assert close(value, breakpoint_inner(seq_y, seq_x))
    assert close(value_norm, breakpoint_inner(seq_y, seq_x, normal=True))
    assert breakpoint_inner(seq_x, list(seq_x)) == 0.0
    assert breakpoint_inner(seq_x, list(seq_x), normal=True) == 0.0
    assert (
        breakpoint_inner(seq_x, seq_z) <= value + breakpoint_inner(seq_y, seq_z) + 1e-9
    )
    assert close(
        order.breakpoint_simil(seq_x, seq_y, boundaries=False), 1.0 - value_norm
    )


# Permutations of distinct items, for the metric claims that hold for them
permutations = st.lists(st.integers(0, 7), unique=True, max_size=8)


@given(perm=permutations, data=st.data())
@settings(max_examples=300, deadline=None)
def test_breakpoint_without_boundaries_permutations(perm, data):
    # For permutations of the same items (at least two), only identical
    # orders score zero
    other = data.draw(st.permutations(perm))
    if len(perm) >= 2 and list(other) != list(perm):
        assert breakpoint_inner(perm, other) > 0.0


@given(seq_x=sequences, seq_y=sequences)
@settings(max_examples=300, deadline=None)
def test_kendall_tau_simil(seq_x, seq_y):
    value = order.kendall_tau_simil(seq_x, seq_y)
    value_norm = order.kendall_tau_simil(seq_x, seq_y, normal=True)
    shared, _ = order.restrict_to_shared(seq_x, seq_y)

    if len(shared) < 2:
        assert math.isnan(value) and math.isnan(value_norm)
    else:
        assert isinstance(value, float)
        assert -1.0 <= value <= 1.0
        assert 0.0 <= value_norm <= 1.0
        assert close(value, order.kendall_tau_simil(seq_y, seq_x))
        assert close(value_norm, (1.0 + value) / 2.0)


@given(seq=sequences)
@settings(max_examples=200, deadline=None)
def test_kendall_tau_simil_identical(seq):
    if len(seq) >= 2:
        assert order.kendall_tau_simil(seq, list(seq)) == 1.0


@given(perm=st.lists(st.integers(0, 20), unique=True, min_size=2, max_size=12))
@settings(max_examples=200, deadline=None)
def test_kendall_tau_simil_reversed(perm):
    assert order.kendall_tau_simil(perm, perm[::-1]) == -1.0
