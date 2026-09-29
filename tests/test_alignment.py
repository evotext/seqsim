"""
test_alignment
==============

Tests for the `alignment` module of the `seqsim` package.

The dynamic programming implementations are compared against brute-force
enumerations of all possible alignments of short sequences.
"""

# Import Python standard libraries
import pytest
from hypothesis import given, settings, strategies as st

# Import the library being tested
from seqsim import alignment, edit

short = st.text(alphabet="abc", max_size=4)


def _alignments(len_x, len_y):
    """
    Yields all alignments as sequences of operations: "M" (aligned pair), "D"
    (element of x against a gap), and "I" (element of y against a gap).
    """

    if len_x == 0 and len_y == 0:
        yield ()
        return
    if len_x and len_y:
        for rest in _alignments(len_x - 1, len_y - 1):
            yield rest + ("M",)
    if len_x:
        for rest in _alignments(len_x - 1, len_y):
            yield rest + ("D",)
    if len_y:
        for rest in _alignments(len_x, len_y - 1):
            yield rest + ("I",)


def _alignment_value(seq_x, seq_y, ops, pair_value, gap_open, gap_extend):
    """
    Returns the sum of pair values and the sum of (positive) gap penalties.
    """

    i = j = 0
    pairs = gaps = 0.0
    previous = None
    for op in ops:
        if op == "M":
            pairs += pair_value(seq_x[i], seq_y[j])
            i += 1
            j += 1
        else:
            gaps += gap_extend + (gap_open if op != previous else 0.0)
            if op == "D":
                i += 1
            else:
                j += 1
        previous = op

    return pairs, gaps


def _brute_nw(seq_x, seq_y, gap_open, gap_extend):
    costs = []
    for ops in _alignments(len(seq_x), len(seq_y)):
        pairs, gaps = _alignment_value(
            seq_x, seq_y, ops, alignment._unit_cost, gap_open, gap_extend
        )
        costs.append(pairs + gaps)
    return min(costs)


def _brute_sw(seq_x, seq_y, gap_open, gap_extend):
    best = 0.0
    for a in range(len(seq_x)):
        for b in range(a + 1, len(seq_x) + 1):
            for c in range(len(seq_y)):
                for d in range(c + 1, len(seq_y) + 1):
                    sub_x, sub_y = seq_x[a:b], seq_y[c:d]
                    for ops in _alignments(len(sub_x), len(sub_y)):
                        pairs, gaps = _alignment_value(
                            sub_x,
                            sub_y,
                            ops,
                            alignment._unit_score,
                            gap_open,
                            gap_extend,
                        )
                        best = max(best, pairs - gaps)
    return best


@given(seq_x=short, seq_y=short, gap_open=st.sampled_from([0.0, 0.5, 2.0]))
@settings(max_examples=150, deadline=None)
def test_nw_matches_brute_force(seq_x, seq_y, gap_open):
    assert alignment.nw_dissim(seq_x, seq_y, gap_open=gap_open) == pytest.approx(
        _brute_nw(seq_x, seq_y, gap_open, 1.0)
    )


@given(seq_x=short, seq_y=short, gap_open=st.sampled_from([0.0, 0.5, 2.0]))
@settings(max_examples=60, deadline=None)
def test_sw_matches_brute_force(seq_x, seq_y, gap_open):
    assert alignment.sw_simil(seq_x, seq_y, gap_open=gap_open) == pytest.approx(
        _brute_sw(seq_x, seq_y, gap_open, 1.0)
    )


@given(
    seq_x=st.text(alphabet="abc", max_size=8), seq_y=st.text(alphabet="abc", max_size=8)
)
@settings(max_examples=150, deadline=None)
def test_nw_defaults_are_levenshtein(seq_x, seq_y):
    assert alignment.nw_dissim(seq_x, seq_y) == edit.levenshtein_dist(seq_x, seq_y)


def test_nw_custom_costs():
    # Treat vowels as interchangeable at a small cost
    def cost(a, b):
        if a == b:
            return 0.0
        if a in "aeiou" and b in "aeiou":
            return 0.1
        return 1.0

    assert alignment.nw_dissim("kitten", "kittan", sub_cost=cost) == pytest.approx(0.1)
    assert alignment.nw_dissim("kitten", "kittan") == 1.0


def test_nw_normal():
    assert alignment.nw_dissim("abc", "xyz", normal=True) == 0.5
    assert alignment.nw_dissim("", "", normal=True) == 0.0
    assert alignment.nw_dissim("abc", "", normal=True) == 1.0


def test_nw_invalid():
    with pytest.raises(ValueError):
        alignment.nw_dissim("abc", "abd", gap_extend=0.0)
    with pytest.raises(ValueError):
        alignment.nw_dissim("abc", "abd", gap_open=-1.0)
    with pytest.raises(ValueError):
        alignment.nw_dissim("abc", "abd", sub_cost=lambda a, b: -1.0)


def test_sw_local():
    assert alignment.sw_simil("XXXabcdYYY", "ZZabcdZZ") == 4.0
    assert alignment.sw_simil("XXXabcdYYY", "ZZabcdZZ", normal=True) == 0.4
    assert alignment.sw_simil("abc", "abc", normal=True) == 1.0
    assert alignment.sw_simil("", "", normal=True) == 1.0
    assert alignment.sw_simil("abc", "", normal=True) == 0.0


def test_monge_elkan():
    titles_x = ["vita antonii", "passio sanctae agnetis"]
    titles_y = ["passio s. agnetis", "vita antonij", "sermo de nativitate"]

    simil = alignment.monge_elkan_simil(titles_x, titles_y)
    assert simil == alignment.monge_elkan_simil(titles_y, titles_x)
    assert 0.0 < simil < 1.0
    assert alignment.monge_elkan_simil(titles_x, list(titles_x)) == 1.0
    assert alignment.monge_elkan_simil([], []) == 1.0
    assert alignment.monge_elkan_simil(titles_x, []) == 0.0

    # Custom inner similarity: exact match only
    exact = alignment.monge_elkan_simil(
        titles_x, titles_y, inner=lambda a, b: float(a == b)
    )
    assert exact == 0.0
