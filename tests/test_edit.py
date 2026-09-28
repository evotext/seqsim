"""
test_edit
=========

Tests for the `edit` module of the `seqsim` package.

Expected values were verified by hand or against independent references (see
the comments of each test). Generic properties (symmetry, identity, triangle
inequality, ranges, and empty sequences) are tested for all methods in
`test_properties.py`.
"""

# Import Python standard libraries
import pytest

# Import the library being tested
from seqsim import edit

# The same four pairs are used for most methods
PAIRS = [
    ("kitten", "sitting"),
    ((1, 2, 3), [1, 2, 3]),
    ((1, 2, 3, 4, 5), (1, 2, 4, 3, 6, 7)),
    ((1, 2, 3), ["a", "b", "c", "d"]),
]


def _cases(values):
    return [(*pair, *value) for pair, value in zip(PAIRS, values)]


# Classic textbook values (kitten/sitting = 3), identical, one transposition
# plus two edits, and fully disjoint sequences
@pytest.mark.parametrize(
    "seq_x,seq_y,expected,expected_norm",
    _cases([(3.0, 3 / 7), (0.0, 0.0), (3.0, 0.5), (4.0, 1.0)]),
)
def test_levenshtein(seq_x, seq_y, expected, expected_norm):
    assert edit.levenshtein_dist(seq_x, seq_y) == expected
    assert edit.levenshtein_dist(seq_x, seq_y, normal=True) == pytest.approx(
        expected_norm
    )


@pytest.mark.parametrize(
    "seq_x,seq_y,expected",
    [
        ("kitten", "sitting", 3.0),
        ("ca", "abc", 2.0),  # transposition plus insertion; OSA gives 3
        ("ab", "ba", 1.0),
        ("abcdef", "badcfe", 3.0),
        ("a cat", "an act", 2.0),
        ((1, 2, 3, 4, 5), (1, 2, 4, 3, 6, 7), 3.0),
        ([], [1, 2], 2.0),
    ],
)
def test_damerau(seq_x, seq_y, expected):
    assert edit.damerau_dist(seq_x, seq_y) == expected
    assert edit.damerau_dist(seq_y, seq_x) == expected


@pytest.mark.parametrize(
    "seq_x,seq_y,expected",
    [
        ("kitten", "sitting", 3.0),
        ("ca", "abc", 3.0),  # the restriction forbids editing "ac" again
        ("ab", "ba", 1.0),
        ("abcdef", "badcfe", 3.0),
        ((1, 2, 3, 4, 5), (1, 2, 4, 3, 6, 7), 3.0),
    ],
)
def test_osa(seq_x, seq_y, expected):
    assert edit.osa_dissim(seq_x, seq_y) == expected
    assert edit.osa_dissim(seq_y, seq_x) == expected


def test_damerau_brute_force():
    """
    Compare the Damerau-Levenshtein distance against a breadth-first search.
    """

    import itertools

    def neighbours(seq, alphabet):
        for i in range(len(seq) + 1):
            for c in alphabet:
                yield seq[:i] + c + seq[i:]
        for i in range(len(seq)):
            yield seq[:i] + seq[i + 1 :]
            for c in alphabet:
                yield seq[:i] + c + seq[i + 1 :]
        for i in range(len(seq) - 1):
            yield seq[:i] + seq[i + 1] + seq[i] + seq[i + 2 :]

    def bfs(source, target, alphabet):
        frontier, seen, dist = {source}, {source}, 0
        while target not in frontier:
            dist += 1
            frontier = {
                n
                for s in frontier
                for n in neighbours(s, alphabet)
                if n not in seen and len(n) <= 5
            }
            seen |= frontier
        return dist

    alphabet = "abc"
    words = [
        "".join(w) for n in range(4) for w in itertools.product(alphabet, repeat=n)
    ]
    for seq_x in words[::3]:
        for seq_y in words[::5]:
            assert edit.damerau_dist(seq_x, seq_y) == bfs(seq_x, seq_y, alphabet)


@pytest.mark.parametrize(
    "seq_x,seq_y,max_del_len,expected",
    [
        ("kitten", "sitting", 5, 3.0),
        ("abcdeXXXXXfghij", "abcdefghij", 5, 1.0),
        ("abcdeXXXXXXfghij", "abcdefghij", 5, 2.0),
        ("abcdeXXXXXfghij", "abcdefghij", 1, 5.0),
        ("abcXdef", "abcdef", 1, 1.0),
        ("abcdefgh", "a", 5, 2.0),  # symmetric block insertion/deletion
        ("a", "abcdefgh", 5, 2.0),
        ((1, 2, 3), ["a", "b", "c", "d"], 5, 2.0),  # delete + insert blocks
    ],
)
def test_bulk_delete(seq_x, seq_y, max_del_len, expected):
    assert edit.bulk_delete_dist(seq_x, seq_y, max_del_len=max_del_len) == expected


def test_bulk_delete_is_levenshtein_with_unit_blocks():
    for seq_x, seq_y in PAIRS:
        assert edit.bulk_delete_dist(
            seq_x, seq_y, max_del_len=1
        ) == edit.levenshtein_dist(seq_x, seq_y)


@pytest.mark.parametrize(
    "seq_x,seq_y,expected",
    [
        # Fragile positions of "kitten" (m=6): 1 and 6; of "sitting" (n=7):
        # 1 and 7. Two substitutions plus inserting the final "g" at half cost.
        ("kitten", "sitting", 2.5),
        # Deleting the first element (fragile) costs half
        ("abcdefghij", "bcdefghij", 0.5),
        ("bcdefghij", "abcdefghij", 0.5),
        # Deleting an interior element costs one
        ("abcdefghij", "abcdfghij", 1.0),
    ],
)
def test_fragile_ends(seq_x, seq_y, expected):
    assert edit.fragile_ends_dissim(seq_x, seq_y) == expected


def test_fragile_ends_without_fragile_regions():
    for seq_x, seq_y in PAIRS:
        assert edit.fragile_ends_dissim(
            seq_x, seq_y, frag_start=0.0, frag_end=0.0
        ) == edit.levenshtein_dist(seq_x, seq_y)


@pytest.mark.parametrize(
    "seq_x,seq_y,kwargs,expected",
    [
        ("kitten", "sitting", {}, 2.5),
        ("abcdeXXXXXfghij", "abcdefghij", {}, 1.0),
        # A block entirely in the fragile start region costs half
        ("XXabcdefghijklmnopqr", "abcdefghijklmnopqr", {"frag_start": 20.0}, 0.5),
        # With no fragile regions this is the "bulk delete" distance
        ("abcdeXXXXXfghij", "abcdefghij", {"frag_start": 0.0, "frag_end": 0.0}, 1.0),
        (
            "abcXdef",
            "abcdef",
            {"frag_start": 0.0, "frag_end": 0.0, "max_del_len": 1},
            1.0,
        ),
    ],
)
def test_stemmatological(seq_x, seq_y, kwargs, expected):
    assert edit.stemmatological_dissim(seq_x, seq_y, **kwargs) == expected


@pytest.mark.parametrize("max_del_len", [0, -1, 1.5, True, None])
def test_invalid_max_del_len(max_del_len):
    with pytest.raises(ValueError):
        edit.bulk_delete_dist("abc", "abd", max_del_len=max_del_len)
    with pytest.raises(ValueError):
        edit.stemmatological_dissim("abc", "abd", max_del_len=max_del_len)


@pytest.mark.parametrize("frag", [-1.0, 100.1])
def test_invalid_frag(frag):
    with pytest.raises(ValueError):
        edit.fragile_ends_dissim("abc", "abd", frag_start=frag)
    with pytest.raises(ValueError):
        edit.stemmatological_dissim("abc", "abd", frag_end=frag)


# Jaro and Jaro-Winkler values match the reference implementation in
# `textdistance` 4.5 and the textbook examples (e.g., "MARTHA"/"MARHTA")
@pytest.mark.parametrize(
    "seq_x,seq_y,jaro,jaro_winkler",
    [
        ("MARTHA", "MARHTA", 1 - 0.944444, 1 - 0.961111),
        ("DIXON", "DICKSONX", 1 - 0.766667, 1 - 0.813333),
        ("kitten", "sitting", 0.253968, 0.253968),
        ((1, 2, 3), [1, 2, 3], 0.0, 0.0),
        ((1, 2, 3, 4, 5), (1, 2, 4, 3, 6, 7), 0.261111, 0.208889),
        ((1, 2, 3), ["a", "b", "c", "d"], 1.0, 1.0),
    ],
)
def test_jaro(seq_x, seq_y, jaro, jaro_winkler):
    assert edit.jaro_dissim(seq_x, seq_y) == pytest.approx(jaro, abs=1e-6)
    assert edit.jaro_winkler_dissim(seq_x, seq_y) == pytest.approx(
        jaro_winkler, abs=1e-6
    )


# MMCWPA values computed by hand: 1 - sqrt(SSNC) / (len(x) + len(y))
@pytest.mark.parametrize(
    "seq_x,seq_y,ssnc",
    [
        ("kitten", "sitting", 6**2 + 2**2),  # "itt", then "n"
        ((1, 2, 3), [1, 2, 3], 6**2),
        ("QabRcd", "abcd", 4**2 + 4**2),  # "ab", then "cd"
        ((1, 2, 3), ["a", "b", "c", "d"], 0),
        ("abc", "bcde", 4**2),
    ],
)
def test_mmcwpa(seq_x, seq_y, ssnc):
    expected = 1 - (ssnc**0.5) / (len(seq_x) + len(seq_y))
    assert edit.mmcwpa_dissim(seq_x, seq_y) == pytest.approx(expected)
    assert edit.mmcwpa_dissim(seq_y, seq_x) == pytest.approx(expected)


# Birnbaum scores computed by hand: sum of n * (n + 1) / 2 over the sizes of
# the matching blocks
@pytest.mark.parametrize(
    "seq_x,seq_y,expected,expected_norm",
    [
        ("kitten", "sitting", 7.0, 7 / 28),  # "itt" + "n"
        ((1, 2, 3), [1, 2, 3], 6.0, 1.0),
        ((1, 2, 3, 4, 5), (1, 2, 4, 3, 6, 7), 4.0, 4 / 21),  # (1, 2) + (3,)
        ((1, 2, 3), ["a", "b", "c", "d"], 0.0, 0.0),
        ("abc", "bcde", 3.0, 3 / 10),
    ],
)
def test_birnbaum(seq_x, seq_y, expected, expected_norm):
    assert edit.birnbaum_simil(seq_x, seq_y) == expected
    assert edit.birnbaum_simil(seq_x, seq_y, normal=True) == pytest.approx(
        expected_norm
    )
    assert edit.birnbaum_dissim(seq_x, seq_y) == pytest.approx(1 - expected_norm)


def test_birnbaum_containment_is_not_identity():
    # In 0.3.1, a sequence contained in another had a distance of zero
    assert edit.birnbaum_dissim("ab", "abXXXXXX") > 0.0


def test_birnbaum_simil_returns_float():
    assert isinstance(edit.birnbaum_simil("abc", "abc"), float)
    assert edit.birnbaum_simil("abc", "abc", normal=True) == 1.0
    assert edit.birnbaum_simil((1, 2, 3), [1, 2, 3], normal=True) == 1.0
