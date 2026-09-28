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


@pytest.mark.parametrize("boundaries", [True, False])
def test_iebp_unbiased_for_few_transpositions(boundaries):
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
            estimates.append(
                order.iebp_estimate(list(range(size)), seq, boundaries=boundaries)
            )
        assert abs(sum(estimates) / len(estimates) - k) < 0.5


def test_iebp_without_boundaries():
    # A transposition of interior blocks breaks three interior adjacencies
    assert order.iebp_estimate("abcdefghij", "abcfghdeij", boundaries=False) == 1.0
    # Moving the last item to the start breaks a single interior adjacency,
    # fewer than the 3 (n - 2) / n = 2.4 expected after one random
    # transposition, so the nearest estimate is zero
    assert order.iebp_estimate("abcdefghij", "jabcdefghi", boundaries=False) == 0.0
    assert order.iebp_estimate("ab", "ba", boundaries=False) == 1.0
    assert order._iebp(10, 0, boundaries=False) == 0


@pytest.mark.parametrize(
    "seq_x,seq_y,kwargs,expected",
    [
        ("abcde", "xdcba", {}, (list("abcd"), list("dcba"))),
        ("abc", "xyz", {}, ([], [])),
        ("", "abc", {}, ([], [])),
        # Occurrences: the second "a" of x has no counterpart in y
        ("abab", "ab", {}, (list("ab"), list("ab"))),
        ("abab", "bab", {}, (list("abb"), list("bab"))),
        ("abab", "bab", {"repeats": "first"}, (list("ab"), list("ba"))),
        ("aab", "ba", {"repeats": "first"}, (list("ab"), list("ba"))),
        ([1, (2, 3), None], [None, 1], {}, ([1, None], [None, 1])),
    ],
)
def test_restrict_to_shared(seq_x, seq_y, kwargs, expected):
    assert order.restrict_to_shared(seq_x, seq_y, **kwargs) == expected


def test_restrict_to_shared_reference():
    # Without repeats, both options give the plain restriction
    import random

    rng = random.Random(3)
    for _ in range(200):
        x = rng.sample(range(30), rng.randint(0, 30))
        y = rng.sample(range(30), rng.randint(0, 30))
        shared = set(x) & set(y)
        expected = ([s for s in x if s in shared], [s for s in y if s in shared])
        assert order.restrict_to_shared(x, y) == expected
        assert order.restrict_to_shared(x, y, repeats="first") == expected

    with pytest.raises(ValueError):
        order.restrict_to_shared("ab", "ab", repeats="last")


@pytest.mark.parametrize(
    "seq_x,seq_y,kwargs,expected",
    [
        # Items only in one sequence cost 1 each, as for top-k lists, so that
        # two permutations reversed reach the maximum
        ("abcd", "dcba", {}, 1.0),
        ("ab", "ba", {}, 1.0),
        ("abcd", "abcd", {}, 0.0),
        # One discordant pair out of C(4, 2)
        ("abcd", "abdc", {}, 1 / 6),
        # Disjoint: 1 (a, x), and the pairs with the boundary cost 1 each
        ("a", "x", {}, 1.0),
        # Maximum for these contents: x-only items first, shared reversed:
        # "Xab" vs "ba" gives the maximum, 4 of 4 changeable pairs
        ("Xab", "ba", {}, 1.0),
        ("abX", "ab", {}, 0.25),
        # "b" and "c" absent from y: pair (b, c) costs p, bound counts p
        ("abc", "a", {"p": 0.5}, 2.5 / 4.5),
        ("abc", "a", {"p": 0.0}, 2.0 / 4.0),
        ("", "", {}, 0.0),
        ("", "a", {}, 1.0),
        ("a", "a", {}, 0.0),
    ],
)
def test_kendall_tau_normal(seq_x, seq_y, kwargs, expected):
    value = order.kendall_tau_dissim(seq_x, seq_y, normal=True, **kwargs)
    assert value == pytest.approx(expected)
    assert order.kendall_tau_dissim(seq_y, seq_x, normal=True, **kwargs) == value


def test_kendall_tau_normal_is_attained():
    # The normalizer is the maximum over all orders of the same contents
    alphabet = "abcd"
    seqs = [
        perm
        for n in range(len(alphabet) + 1)
        for comb in itertools.combinations(alphabet, n)
        for perm in itertools.permutations(comb)
    ]
    for p in (0.0, 0.5, 1.0):
        maxima = {}
        for x in seqs:
            for y in seqs:
                key = (frozenset(x), frozenset(y))
                value = order.kendall_tau_dissim(x, y, p=p, normal=True)
                assert 0.0 <= value <= 1.0
                maxima[key] = max(maxima.get(key, 0.0), value)
        for (set_x, set_y), value in maxima.items():
            # Only a single shared item (or nothing) cannot be reordered
            if (set_x, set_y) != (set_x, set_x) or len(set_x) > 1:
                assert value == pytest.approx(1.0)


def test_kendall_tau_normal_permutations():
    # For permutations, the classic normalized Kendall distance
    for target in PERMS[::5]:
        dist = order.kendall_tau_dissim("abcde", target)
        assert order.kendall_tau_dissim("abcde", target, normal=True) == dist / 10


@pytest.mark.parametrize(
    "seq_x,seq_y,expected",
    [
        ("abcd", "abcd", 1.0),
        ("abcd", "dcba", -1.0),
        # One discordant pair of six: 1 - 4 / 12
        ("abcd", "abdc", 2 / 3),
        # Only shared items count: "X" and "Y" are ignored
        ("aXbcd", "abdcY", 2 / 3),
        # Two discordant pairs of three: 1 - 8 / 6
        ("abc", "cab", -1 / 3),
        ("ab", "ba", -1.0),
    ],
)
def test_kendall_tau_simil(seq_x, seq_y, expected):
    assert order.kendall_tau_simil(seq_x, seq_y) == pytest.approx(expected)
    assert order.kendall_tau_simil(seq_y, seq_x) == pytest.approx(expected)
    assert order.kendall_tau_simil(seq_x, seq_y, normal=True) == pytest.approx(
        (1 + expected) / 2
    )


def test_kendall_tau_simil_undefined():
    import math

    for seq_x, seq_y in [("", ""), ("a", "a"), ("abc", "xyz"), ("ab", "b")]:
        assert math.isnan(order.kendall_tau_simil(seq_x, seq_y))
        assert math.isnan(order.kendall_tau_simil(seq_x, seq_y, normal=True))


def test_kendall_tau_simil_matches_distance():
    # tau = 1 - 4 d / (n (n - 1)), d the Kendall distance of the reduced orders
    import random

    rng = random.Random(7)
    for _ in range(200):
        x = rng.sample(range(40), rng.randint(0, 40))
        y = rng.sample(range(40), rng.randint(0, 40))
        rx, ry = order.restrict_to_shared(x, y)
        n = len(rx)
        if n < 2:
            continue
        d = order.kendall_tau_dissim(rx, ry)
        assert order.kendall_tau_simil(x, y) == pytest.approx(1 - 4 * d / (n * (n - 1)))


def test_count_inversions():
    for perm in itertools.permutations(range(6)):
        expected = sum(
            1 for i, j in itertools.combinations(range(6), 2) if perm[i] > perm[j]
        )
        assert order._count_inversions(list(perm)) == expected


@pytest.mark.parametrize(
    "seq_x,seq_y,expected",
    [
        # a/b, b/c kept; c/d, d/e, e/f broken
        ("abcdef", "abcfed", 3.0),
        ("abcdef", "abcdef", 0.0),
        # Moving an end item costs one, not two
        ("abcdef", "fabcde", 1.0),
        ("abc", "xyz", 2.0),
        # No adjacencies without boundaries
        ("", "", 0.0),
        ("a", "b", 0.0),
        ("", "abc", 1.0),
    ],
)
def test_breakpoint_without_boundaries(seq_x, seq_y, expected):
    assert order.breakpoint_dissim(seq_x, seq_y, boundaries=False) == expected
    assert order.breakpoint_dissim(seq_y, seq_x, boundaries=False) == expected


def test_breakpoint_without_boundaries_normal():
    assert (
        order.breakpoint_dissim("abcdef", "abcfed", boundaries=False, normal=True)
        == 0.6
    )
    assert order.breakpoint_dissim("abc", "xyz", boundaries=False, normal=True) == 1.0
    assert order.breakpoint_dissim("a", "b", boundaries=False, normal=True) == 0.0
    assert order.breakpoint_dissim("", "abc", boundaries=False, normal=True) == 1.0


def test_breakpoint_simil():
    assert order.breakpoint_simil("abcdef", "abcfed", boundaries=False) == 0.4
    assert order.breakpoint_simil("abcdef", "abcfed") == 3 / 7
    assert order.breakpoint_simil("", "") == 1.0
    assert order.breakpoint_simil("a", "b", boundaries=False) == 1.0
    assert order.breakpoint_simil("abc", "xyz") == 0.0


def test_breakpoint_simil_is_share_of_adjacencies():
    # The share of adjacencies of the shared items preserved, exactly
    import random

    rng = random.Random(11)
    for _ in range(500):
        pool = range(rng.randint(2, 60))
        x = rng.sample(pool, rng.randint(2, len(pool)))
        y = rng.sample(pool, rng.randint(2, len(pool)))
        rx, ry = order.restrict_to_shared(x, y)
        if len(rx) < 2:
            continue
        shared = set(zip(rx, rx[1:])) & set(zip(ry, ry[1:]))
        expected = len(shared) / (len(rx) - 1)
        assert order.breakpoint_simil(rx, ry, boundaries=False) == expected
        assert 1 - order.breakpoint_dissim(
            rx, ry, boundaries=False, normal=True
        ) == pytest.approx(expected)
