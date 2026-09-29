"""
test_properties
===============

Property-based tests (with `hypothesis`) of the claims that each measure
declares in its registration (see `seqsim.measures()`):

  * all measures: float results, and normalized values in range [0..1]
    (except for estimators, whose normalized values are rates);
  * `normal`: no effect where declared so;
  * symmetry, where declared;
  * identical sequences: `0.0` for distances, dissimilarities and
    estimators, the largest value for similarities;
  * identity of indiscernibles: tested where claimed, and the stored
    example verified where it fails;
  * triangle inequality: tested where claimed (under its condition, where
    conditional), and the stored counterexample verified where it fails;
  * empty sequences, following the declared empty rule.

A measure added with its declaration is thus tested without changes here.
A few tests at the end cover options and relations between measures that
the declarations do not describe.
"""

# Import Python standard libraries
import math

import pytest
from hypothesis import given, settings, strategies as st

# Import the library being tested
import seqsim
from seqsim import edit, order

INFOS = seqsim.measures()

# Small alphabets make shared elements, and thus interesting cases, likely;
# elements of different types are mixed to exercise arbitrary hashables
elements = st.sampled_from(["a", "b", "c", 1, (2, 3), None])
seqs = st.lists(elements, max_size=8)
strings = st.text(alphabet="abc", max_size=8)
sequences = st.one_of(seqs, strings, seqs.map(tuple))

# Measures whose elements must themselves be sequences
STRATEGIES = {
    "alignment.monge_elkan_simil": st.lists(strings, max_size=4),
}

# Options under which conditional triangle inequalities hold
CONDITIONS = {
    "alignment.nw_dissim": {"gap_open": 0.0},
    "order.footrule_dissim": {"ell": 10},
}

# Measures without a value for some pairs (fewer than two shared elements)
UNDEFINED = {"order.kendall_tau_simil"}


def select(predicate):
    return [pytest.param(info, id=info.name) for info in INFOS if predicate(info)]


def inputs(info):
    return STRATEGIES.get(info.name, sequences)


def close(value, target):
    return math.isclose(value, target, abs_tol=1e-9)


def largest(info):
    """The value of identical sequences: 0.0, or 1.0 for similarities."""

    return 1.0 if info.kind == "simil" else 0.0


def values(info, seq_x, seq_y, **kwargs):
    """The raw and, if accepted, the normalized value."""

    result = [info.func(seq_x, seq_y, **kwargs)]
    if info.accepts_normal:
        result.append(info.func(seq_x, seq_y, normal=True, **kwargs))
    return result


def test_declarations():
    names = [info.name for info in INFOS]
    keys = [info.key for info in INFOS if info.key]
    assert len(names) == len(set(names))
    assert len(keys) == len(set(keys))
    assert keys and sorted(keys) == list(seqsim.METHODS)
    assert set(CONDITIONS) == {
        info.name for info in INFOS if info.triangle == "conditional"
    }
    for info in INFOS:
        suffix = {"dist": "_dist", "dissim": "_dissim", "simil": "_simil"}
        if info.kind in suffix:
            assert info.name.endswith(suffix[info.kind]), info.name
        assert info.func.info is info
        # Measures available in `distance()` are dissimilarities
        assert info.key is None or info.kind in ("dist", "dissim")


@given(data=st.data())
@settings(max_examples=150, deadline=None)
@pytest.mark.parametrize("info", select(lambda info: True))
def test_range_and_type(info, data):
    seq_x, seq_y = data.draw(inputs(info)), data.draw(inputs(info))
    value, *normalized = values(info, seq_x, seq_y)

    assert isinstance(value, float)
    if info.name in UNDEFINED and math.isnan(value):
        return
    if info.kind != "simil" or info.raw_range.startswith("0"):
        assert value >= 0.0
    for value_norm in normalized:
        assert isinstance(value_norm, float)
        assert value_norm >= 0.0
        assert info.kind == "estimator" or value_norm <= 1.0


@given(data=st.data())
@settings(max_examples=100, deadline=None)
@pytest.mark.parametrize(
    "info", select(lambda info: info.accepts_normal and not info.normal_effect)
)
def test_normal_without_effect(info, data):
    seq_x, seq_y = data.draw(inputs(info)), data.draw(inputs(info))
    assert info.func(seq_x, seq_y) == info.func(seq_x, seq_y, normal=True)


@given(data=st.data())
@settings(max_examples=150, deadline=None)
@pytest.mark.parametrize("info", select(lambda info: info.symmetric))
def test_symmetry(info, data):
    seq_x, seq_y = data.draw(inputs(info)), data.draw(inputs(info))
    for value_xy, value_yx in zip(
        values(info, seq_x, seq_y), values(info, seq_y, seq_x)
    ):
        assert close(value_xy, value_yx) or (
            math.isnan(value_xy) and math.isnan(value_yx)
        )


@given(data=st.data())
@settings(max_examples=100, deadline=None)
@pytest.mark.parametrize("info", select(lambda info: info.identical_zero))
def test_identical(info, data):
    seq = data.draw(inputs(info))
    value, *normalized = values(info, seq, list(seq))

    if info.name in UNDEFINED and math.isnan(value):
        return
    if info.kind != "simil":
        assert value == 0.0
    for value_norm in normalized:
        assert value_norm == largest(info)
    if not info.accepts_normal:
        assert value == largest(info)


@given(data=st.data())
@settings(max_examples=150, deadline=None)
@pytest.mark.parametrize("info", select(lambda info: info.identity == "yes"))
def test_identity_of_indiscernibles(info, data):
    seq_x, seq_y = data.draw(inputs(info)), data.draw(inputs(info))
    if list(seq_x) != list(seq_y):
        for value in values(info, seq_x, seq_y):
            assert value > 0.0


@pytest.mark.parametrize("info", select(lambda info: info.identity == "no"))
def test_identity_counterexample(info):
    seq_x, seq_y = info.identity_example
    assert list(seq_x) != list(seq_y)
    assert info.func(seq_x, seq_y) == 0.0


@given(data=st.data())
@settings(max_examples=300, deadline=None)
@pytest.mark.parametrize(
    "info", select(lambda info: info.triangle in ("yes", "conditional"))
)
def test_triangle_inequality(info, data):
    seq_x, seq_y, seq_z = (data.draw(inputs(info)) for _ in range(3))
    kwargs = CONDITIONS.get(info.name, {})

    def dist(seq_a, seq_b):
        return info.func(seq_a, seq_b, **kwargs)

    assert dist(seq_x, seq_z) <= dist(seq_x, seq_y) + dist(seq_y, seq_z) + 1e-9


@pytest.mark.parametrize("info", select(lambda info: info.triangle == "no"))
def test_triangle_counterexample(info):
    seq_x, seq_y, seq_z = info.triangle_example
    dist = info.func
    assert dist(seq_x, seq_z) > dist(seq_x, seq_y) + dist(seq_y, seq_z) + 1e-9


@given(data=st.data())
@settings(max_examples=50, deadline=None)
@pytest.mark.parametrize(
    "info",
    select(
        lambda info: info.kind in ("dist", "dissim", "simil")
        and info.name not in UNDEFINED
    ),
)
def test_empty_sequences(info, data):
    seq = data.draw(inputs(info).filter(len))
    empty = type(seq)()
    both_empty = values(info, empty, empty)
    one_empty = values(info, seq, empty) + values(info, empty, seq)

    if info.kind == "simil":
        # The raw value of a similarity may be a score, zero for no elements
        assert both_empty[-1] == 1.0
        assert all(value == 0.0 for value in one_empty)
    else:
        assert all(value == 0.0 for value in both_empty)
        assert all(value > 0.0 for value in one_empty)
        if info.empty == "max":
            assert all(value == 1.0 for value in one_empty)


# Options and relations that the declarations do not describe


@pytest.mark.parametrize("max_del_len", [1, 2, 3])
@given(seq_x=strings, seq_y=strings, seq_z=strings)
@settings(max_examples=200, deadline=None)
def test_bulk_delete_triangle_inequality(max_del_len, seq_x, seq_y, seq_z):
    def dist(a, b):
        return edit.bulk_delete_dist(a, b, max_del_len=max_del_len)

    assert dist(seq_x, seq_z) <= dist(seq_x, seq_y) + dist(seq_y, seq_z)


@given(seq_x=sequences, seq_y=sequences)
@settings(max_examples=150, deadline=None)
def test_birnbaum_simil_identity(seq_x, seq_y):
    simil_norm = edit.birnbaum_simil(seq_x, seq_y, normal=True)
    assert (simil_norm == 1.0) == (list(seq_x) == list(seq_y))


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
        assert -1.0 <= value <= 1.0
        assert close(value_norm, (1.0 + value) / 2.0)


@given(perm=st.lists(st.integers(0, 20), unique=True, min_size=2, max_size=12))
@settings(max_examples=200, deadline=None)
def test_kendall_tau_simil_reversed(perm):
    assert order.kendall_tau_simil(perm, perm[::-1]) == -1.0
