"""
test_order
==========

Tests for the `order` module of the `seqsim` package.

For permutations, distances are compared against breadth-first searches over
the operations that define them.
"""

# Import Python standard libraries
import itertools

import pytest

# Import the library being tested
from seqsim import order


def _bfs(source, target, neighbours):
    """
    Returns the minimum number of operations to transform `source` into `target`.
    """

    frontier, seen, dist = {source}, {source}, 0
    while target not in frontier:
        dist += 1
        frontier = {n for s in frontier for n in neighbours(s) if n not in seen}
        seen |= frontier
    return dist


def _moves(perm):
    for i in range(len(perm)):
        rest = perm[:i] + perm[i + 1 :]
        for j in range(len(perm)):
            yield rest[:j] + (perm[i],) + rest[j:]


def _swaps(perm):
    for i, j in itertools.combinations(range(len(perm)), 2):
        new = list(perm)
        new[i], new[j] = new[j], new[i]
        yield tuple(new)


PERMS = list(itertools.permutations("abcde"))


@pytest.mark.parametrize("target", PERMS[::7])
def test_ulam_permutations(target):
    source = tuple("abcde")
    assert order.ulam_dist(source, target) == _bfs(source, target, _moves)


@pytest.mark.parametrize("target", PERMS[::7])
def test_cayley_permutations(target):
    source = tuple("abcde")
    assert order.cayley_dissim(source, target) == _bfs(source, target, _swaps)


@pytest.mark.parametrize(
    "seq_x,seq_y,expected",
    [
        ("abcdef", "bcdefa", 1.0),
        ("abcdef", "bcdXefa", 2.0),  # one move, one insertion
        ("abc", "xyz", 6.0),  # three deletions, three insertions
        ("aab", "aba", 1.0),  # occurrences are matched in order
        ("", "", 0.0),
    ],
)
def test_ulam(seq_x, seq_y, expected):
    assert order.ulam_dist(seq_x, seq_y) == expected


def test_ulam_normal():
    assert order.ulam_dist("abc", "xyz", normal=True) == 1.0
    assert order.ulam_dist("abcd", "dcba", normal=True) == 0.75
    assert order.ulam_dist("", "", normal=True) == 0.0


@pytest.mark.parametrize(
    "seq_x,seq_y,expected",
    [
        # Start/a, a/b, ..., f/end: c/d, d/e, e/f, f/end are broken
        ("abcdef", "abcfed", 4.0),
        ("abcdef", "abcdef", 0.0),
        ("abc", "xyz", 4.0),
        ("abacada", "acabada", 0.0),  # same adjacencies, different sequences
        ("", "a", 1.5),
    ],
)
def test_breakpoint(seq_x, seq_y, expected):
    assert order.breakpoint_dissim(seq_x, seq_y) == expected
    assert order.breakpoint_dissim(seq_y, seq_x) == expected


def test_breakpoint_normal():
    assert order.breakpoint_dissim("abc", "xyz", normal=True) == 1.0
    assert order.breakpoint_dissim("", "", normal=True) == 0.0


def _block_interchanges(perm):
    n = len(perm)
    for i in range(n):
        for j in range(i + 1, n + 1):
            for k in range(j, n):
                for m in range(k + 1, n + 1):
                    yield perm[:i] + perm[k:m] + perm[j:k] + perm[i:j] + perm[m:]


@pytest.mark.parametrize("target", list(itertools.permutations("abcdef"))[::37])
def test_block_interchange_permutations(target):
    source = tuple("abcdef")
    assert order.block_interchange_dissim(source, target) == _bfs(
        source, target, _block_interchanges
    )


def _ulam_neighbours(seq, alphabet="ab"):
    # Moves, insertions, and deletions of single elements
    for i in range(len(seq)):
        rest = seq[:i] + seq[i + 1 :]
        yield rest
        for j in range(len(seq)):
            yield rest[:j] + seq[i] + rest[j:]
    for i in range(len(seq) + 1):
        for c in alphabet:
            yield seq[:i] + c + seq[i:]


def test_ulam_optimal_with_repetitions():
    # With repeated elements, matching occurrences in order is not optimal:
    # "ba" -> "aba" needs a single insertion
    assert order.ulam_dist("ba", "aba") == 1.0
    words = ["".join(w) for n in range(5) for w in itertools.product("ab", repeat=n)]
    for seq_x in words[::3]:
        for seq_y in words[::4]:
            expected = _bfs(
                seq_x, seq_y, lambda s: (n for n in _ulam_neighbours(s) if len(n) <= 6)
            )
            assert order.ulam_dist(seq_x, seq_y) == expected


@pytest.mark.parametrize(
    "seq_x,seq_y,kwargs,expected",
    [
        ("abcd", "bacd", {}, 1.0),
        ("abcd", "dcba", {}, 6.0),
        # Top-k lists (Fagin et al.): "ab" vs "ac" share "a"; (b, c) are
        # each in one list only (1); (b, END) and (c, END) add 1 each
        ("ab", "ac", {}, 3.0),
        # Both "b" and "c" absent from the other list, in the same list
        ("abc", "a", {"p": 0.5}, 2.5),
        ("abc", "a", {"p": 1.0}, 3.0),
        ("", "", {}, 0.0),
        ("", "a", {}, 1.0),
    ],
)
def test_kendall_tau(seq_x, seq_y, kwargs, expected):
    assert order.kendall_tau_dissim(seq_x, seq_y, **kwargs) == expected
    assert order.kendall_tau_dissim(seq_y, seq_x, **kwargs) == expected


def test_kendall_tau_permutations():
    # For permutations, the number of discordant pairs is the minimum number
    # of swaps of adjacent items
    def adjacent_swaps(perm):
        for i in range(len(perm) - 1):
            yield perm[:i] + (perm[i + 1], perm[i]) + perm[i + 2 :]

    source = tuple("abcde")
    for target in PERMS[::5]:
        assert order.kendall_tau_dissim(source, target) == _bfs(
            source, target, adjacent_swaps
        )

    assert order.kendall_tau_dissim("abcd", "dcba", normal=True) == 0.6
    with pytest.raises(ValueError):
        order.kendall_tau_dissim("ab", "ba", p=1.5)


@pytest.mark.parametrize(
    "seq_x,seq_y,kwargs,expected",
    [
        ("abcd", "bacd", {}, 2.0),
        ("abcd", "dcba", {}, 8.0),
        # Missing items at position 3: |2 - 3| + |3 - 2|
        ("ab", "ac", {}, 2.0),
        ("ab", "ac", {"ell": 10}, 16.0),
        ("", "", {}, 0.0),
    ],
)
def test_footrule(seq_x, seq_y, kwargs, expected):
    assert order.footrule_dissim(seq_x, seq_y, **kwargs) == expected
    assert order.footrule_dissim(seq_y, seq_x, **kwargs) == expected


def test_footrule_normal_and_errors():
    assert order.footrule_dissim("abc", "xyz", normal=True) == 1.0
    assert order.footrule_dissim("abc", "abc", normal=True) == 0.0
    with pytest.raises(ValueError):
        order.footrule_dissim("abc", "ab", ell=3)


def test_iebp():
    assert order.iebp_estimate("abcdefghij", "abcdefghij") == 0.0
    assert order.iebp_estimate("abcdefghij", "abcfghdeij") == 1.0
    # Only shared items are considered
    assert order.iebp_estimate("abcdefghij", "abXcdYefghij") == 0.0
    assert order.iebp_estimate("ab", "ba") == 1.0
    assert order.iebp_estimate("", "") == 0.0
    assert order.iebp_estimate("abcdefghij", "abcfghdeij", normal=True) == 0.1


def test_iebp_unbiased_for_few_transpositions():
    import random

    rng = random.Random("seqsim")
    size = 40
    for k in (1, 3, 6):
        estimates = []
        for _ in range(100):
            seq = list(range(size))
            for _ in range(k):
                i, j, m = sorted(rng.sample(range(size + 1), 3))
                seq = seq[:i] + seq[j:m] + seq[i:j] + seq[m:]
            estimates.append(order.iebp_estimate(list(range(size)), seq))
        assert abs(sum(estimates) / len(estimates) - k) < 0.5
