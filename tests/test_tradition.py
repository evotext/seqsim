"""
test_tradition
==============

Tests for the `tradition` subpackage of the `seqsim` package, with expected
values checked by hand. The equivalence with the reference implementation
of the Apophthegmata project is tested in `test_tradition_reference.py`.
"""

# Import Python standard libraries
import os
import pathlib
import subprocess
import sys

import pytest

# Import the library being tested
from seqsim import tradition
from seqsim.tradition import CharacterMatrix

# Frame and labels


def test_unique():
    assert tradition.unique("abacb") == ["a", "b", "c"]
    assert tradition.unique([]) == []
    assert tradition.unique([1, (2,), 1, None]) == [1, (2,), None]


@pytest.mark.parametrize(
    "sequences,expected",
    [
        ([], []),
        (["abc"], ["a", "b", "c"]),
        (["aba"], ["a", "b"]),
        # "c" goes after "b", the last placed item preceding it in the witness
        (["abd", "bcd"], ["a", "b", "c", "d"]),
        # Items before any placed item go at the start
        (["cd", "abc"], ["a", "b", "c", "d"]),
        # The first witness fixes the order of its items
        (["abc", "cba"], ["a", "b", "c"]),
        # "x" goes after "c", its nearest placed predecessor in the witness,
        # although "c" comes after "a" in the frame
        (["abc", "cxa"], ["a", "b", "c", "x"]),
        (["ab", "xy"], ["x", "y", "a", "b"]),
    ],
)
def test_consensus_order(sequences, expected):
    assert tradition.consensus_order(sequences) == expected


def test_consensus_order_depends_on_order():
    assert tradition.consensus_order(["ab", "ba"]) == ["a", "b"]
    assert tradition.consensus_order(["ba", "ab"]) == ["b", "a"]


def test_proportional_labels():
    # Items at positions 0, 0.5, 1, ... 2.5 paragraphs
    events = [(0, "I"), (2, "II")]
    assert tradition.proportional_labels(6, 3, events) == ["I"] * 4 + ["II"] * 2
    # The comparison is <=: an event at 1.5 labels item 3 (position 1.5)
    assert tradition.proportional_labels(4, 2, [(0, "I"), (1.5, "II")]) == [
        "I",
        "I",
        "I",
        "II",
    ]
    # The last event (in the order given) whose position is reached wins
    assert tradition.proportional_labels(2, 2, [(0, "I"), (0, "II")]) == ["II", "II"]
    assert tradition.proportional_labels(2, 2, [(1, "II"), (0, "I")]) == ["I", "I"]
    assert tradition.proportional_labels(3, 3, [(1, "I")]) == [None, "I", "I"]
    assert tradition.proportional_labels(3, 0, [(0, "I")]) == [None] * 3
    assert tradition.proportional_labels(3, 3, []) == [None] * 3
    assert tradition.proportional_labels(0, 3, [(0, "I")]) == []


def test_proportional_positions():
    assert tradition.proportional_positions(10, 4, [0, 1, 4]) == [0, 2, 10]
    # Python's round: halves to even
    assert tradition.proportional_positions(5, 2, [1]) == [2]
    assert tradition.proportional_positions(7, 2, [1]) == [4]
    assert tradition.proportional_positions(10, 0, [1]) == []
    assert tradition.proportional_positions(10, 4, []) == []


def test_monotone_labels():
    reference = {"a": "I", "b": "I", "c": "II", "d": "I", "e": "II"}
    # Labels cannot go back from II to I: "d" loses its reference label
    assert tradition.monotone_labels("abcde", ["I", "II"], reference, [None] * 5) == [
        "I",
        "I",
        "I",
        "I",
        "II",
    ]
    # Without any information, the first label (ties go to the first maximum)
    assert tradition.monotone_labels("xyz", ["I", "II"], {}, [None] * 3) == ["I"] * 3
    # The prior decides items without a reference label
    assert tradition.monotone_labels("xyz", ["I", "II"], {}, ["I", "II", "II"]) == [
        "I",
        "II",
        "II",
    ]
    # A heavy prior outweighs the reference
    prior = ["II", "II", "II"]
    assert tradition.monotone_labels("abc", ["I", "II"], {"a": "I"}, prior) == [
        "I",
        "II",
        "II",
    ]
    assert tradition.monotone_labels(
        "abc", ["I", "II"], {"a": "I"}, prior, prior_weight=2.0
    ) == ["II", "II", "II"]


def test_monotone_labels_empty():
    assert tradition.monotone_labels("", ["I"], {}, []) == []
    assert tradition.monotone_labels("ab", [], {}, [None, None]) == ["", ""]


def test_reference_labels():
    voters = [
        ("abcdef", 6, [(0, "I"), (3, "II")]),
        ("abcdef", 6, [(0, "I"), (2, "II")]),
        ("abdcef", 6, [(0, "I"), (3, "II")]),
        ("xyz", 3, []),  # without events: ignored
    ]
    label_of, support = tradition.reference_labels(voters, rounds=1)
    assert label_of == {"a": "I", "b": "I", "c": "II", "d": "II", "e": "II", "f": "II"}
    assert list(label_of) == list("abcdef")
    assert support["a"] == 1.0
    assert support["c"] == support["d"] == 2 / 3

    label_of, support = tradition.reference_labels(voters)
    assert set(support.values()) == {1.0}

    assert tradition.reference_labels(voters, rounds=0) == ({}, {})
    assert tradition.reference_labels([]) == ({}, {})


def test_reference_labels_ties():
    # A tie goes to the label voted for first, in the order of the voters
    voters = [("ab", 2, [(0, "II")]), ("ab", 2, [(0, "I")])]
    label_of, support = tradition.reference_labels(voters, rounds=1)
    assert label_of == {"a": "II", "b": "II"}
    assert support == {"a": 0.5, "b": 0.5}
    voters = voters[::-1]
    assert tradition.reference_labels(voters, rounds=1)[0] == {"a": "I", "b": "I"}


def test_fill_labels():
    assert tradition.fill_labels("abcde", {"b": "I", "d": "II"}) == {
        "b": "I",
        "d": "II",
        "a": "I",
        "c": "I",
        "e": "II",
    }
    assert tradition.fill_labels("abc", {}) == {}
    # Labels of items outside the frame are kept
    assert tradition.fill_labels("a", {"z": "I"}) == {"z": "I"}


# Coverage

FRAME = list("abcdefghij")
LABEL_OF = dict(zip(FRAME, "IIIIIJJJJJ"))


def test_covered_labels():
    assert tradition.covered_labels("abf", FRAME, LABEL_OF) == {"I"}
    assert tradition.covered_labels("abf", FRAME, LABEL_OF, min_items=1) == {"I", "J"}
    # 2 of 5 items of "I" is 40%
    assert tradition.covered_labels("ab", FRAME, LABEL_OF, min_share=0.4) == {"I"}
    assert tradition.covered_labels("ab", FRAME, LABEL_OF, min_share=0.41) == set()
    # Repeated items and items outside the frame do not count
    assert tradition.covered_labels("aaxy", FRAME, LABEL_OF) == set()
    # Unlabelled items are never covered
    assert tradition.covered_labels("ab", "ab", {}, min_items=1) == set()


def test_coverage():
    cells = tradition.coverage("bcdgh", FRAME, LABEL_OF)
    assert list(cells) == FRAME
    assert cells == {
        "a": None,  # before the first item
        "b": 1,
        "c": 1,
        "d": 1,
        "e": 0,
        "f": 0,
        "g": 1,
        "h": 1,
        "i": None,  # after the last item
        "j": None,
    }


def test_coverage_uncovered_label():
    # Only one item of "J": the label is not covered, its absences are missing
    cells = tradition.coverage("acdh", FRAME, LABEL_OF)
    assert cells["b"] == 0 and cells["e"] == 0
    assert cells["f"] is None and cells["g"] is None
    assert cells["h"] == 1


def test_coverage_gaps():
    # A lacuna before the third item of "bdgi" ("g"): the frame items between
    # "d" and "g" are missing
    cells = tradition.coverage("bdgi", FRAME, LABEL_OF, gaps=[2])
    assert cells["c"] == 0
    assert cells["e"] is None and cells["f"] is None
    assert cells["h"] == 0
    # At the start or the end, a lacuna only affects items already missing
    assert tradition.coverage("bdgi", FRAME, LABEL_OF, gaps=[0, 4]) == (
        tradition.coverage("bdgi", FRAME, LABEL_OF)
    )
    # Gaps index the witness reduced to first occurrences
    assert tradition.coverage("bbdgi", FRAME, LABEL_OF, gaps=[2]) == cells
    # A lacuna between items in reverse frame order loses nothing
    assert tradition.coverage("gdbi", FRAME, LABEL_OF, gaps=[1]) == tradition.coverage(
        "gdbi", FRAME, LABEL_OF
    )


def test_coverage_empty():
    assert tradition.coverage("", FRAME, LABEL_OF) == dict.fromkeys(FRAME)
    assert tradition.coverage("xyz", FRAME, LABEL_OF) == dict.fromkeys(FRAME)
    assert tradition.coverage("ab", [], {}) == {}


def test_merge_coverage():
    part_1 = {"a": 1, "b": 0, "c": None, "d": 0}
    part_2 = {"a": None, "b": None, "c": 1, "d": 1}
    assert tradition.merge_coverage([part_1, part_2], "abcd") == {
        "a": 1,
        "b": 0,
        "c": 1,
        "d": 1,
    }
    assert tradition.merge_coverage([part_1], "abcd") == part_1
    assert tradition.merge_coverage([part_1], "abcd") is not part_1
    assert tradition.merge_coverage([], "ab") == {"a": None, "b": None}


def test_density():
    assert tradition.density({"a": 1, "b": 0, "c": None, "d": 1}) == 2 / 3
    assert tradition.density({"a": None}) == 0.0
    assert tradition.density({}) == 0.0


def test_merge_parts():
    # "B" shares nothing with "A", "C" is mostly contained in it
    parts = [("C", "abx"), ("A", "abcd"), ("B", "efg")]
    merged, separate = tradition.merge_parts(parts)
    assert merged == [("A", "abcd"), ("B", "efg")]
    assert separate == [("C", "abx")]

    # With a high threshold, "C" (2 of 3 items in "A") is merged too
    merged, separate = tradition.merge_parts(parts, max_overlap=0.7)
    assert [label for label, _ in merged] == ["A", "B", "C"]
    assert separate == []

    # Overlap is measured against the merged content so far: "C" overlaps
    # with "B", merged before it
    parts = [("A", "abcdef"), ("B", "ghij"), ("C", "hijk")]
    merged, separate = tradition.merge_parts(parts)
    assert [label for label, _ in merged] == ["A", "B"]
    assert separate == [("C", "hijk")]


def test_merge_parts_ties_and_compatibility():
    # Equal sizes keep the input order: "B" is the largest
    parts = [("B", "abc"), ("A", "def"), ("C", "abg")]
    merged, separate = tradition.merge_parts(parts)
    assert merged == [("A", "def"), ("B", "abc")]
    assert separate == [("C", "abg")]

    # `compatible` compares with the largest part
    def compatible(a, b):
        return {a, b} != {"A", "B"}

    merged, separate = tradition.merge_parts(parts, compatible=compatible)
    assert merged == []
    assert separate == parts


def test_merge_parts_sizes_by_distinct_items():
    # "A" has the most items but the fewest distinct ones: "B" is the largest
    # part, so only parts compatible with it can be merged
    parts = [("A", "aaaaaab"), ("B", "cde"), ("C", "fg")]

    def compatible(a, b):
        return {a, b} == {"A", "C"}

    assert tradition.merge_parts(parts, compatible=compatible) == ([], parts)

    # Without constraints, all parts are merged, in label order
    merged, separate = tradition.merge_parts(parts)
    assert [label for label, _ in merged] == ["A", "B", "C"]
    assert separate == []


def test_merge_parts_nothing_to_merge():
    assert tradition.merge_parts([]) == ([], [])
    assert tradition.merge_parts([("A", "abc")]) == ([], [("A", "abc")])
    parts = [("A", "ab"), ("B", "abc")]
    assert tradition.merge_parts(parts) == ([], parts)
    # Empty parts have no overlap
    assert tradition.merge_parts([("A", "ab"), ("B", "")]) == (
        [("A", "ab"), ("B", "")],
        [],
    )


# Characters


def test_character_matrix():
    matrix = CharacterMatrix(["W1", "W2"], ["a", "b"], [[1, 0], [None, 1]])
    assert matrix.taxa == ("W1", "W2")
    assert matrix.cells == ((1, 0), (None, 1))
    assert matrix.shape == (2, 2)
    assert matrix.column("b") == [0, 1]
    assert matrix.row("W2") == (None, 1)
    assert matrix.as_dict() == {"a": [1, None], "b": [0, 1]}
    assert CharacterMatrix.from_columns(["W1", "W2"], matrix.as_dict()) == matrix
    empty = CharacterMatrix.from_columns(["W1"], {})
    assert empty.shape == (1, 0) and empty.cells == ((),)
    with pytest.raises(ValueError):
        CharacterMatrix(["W1"], ["a"], [])
    with pytest.raises(ValueError):
        CharacterMatrix(["W1"], ["a"], [[1, 0]])
    with pytest.raises(Exception):
        matrix.taxa = ("W3",)


def test_informative():
    assert tradition.informative([1, 1, 0, 0])
    assert not tradition.informative([1, 0, 0, 0, None])
    assert not tradition.informative([None, None])
    assert not tradition.informative([])
    assert tradition.informative([1, 0], min_each=1)
    assert not tradition.informative([1, 1, 0, 0], min_each=3)


def test_content_characters():
    coverages = {
        "W1": {"a": 1, "b": 1, "c": 1, "d": 0},
        "W2": {"a": 1, "b": 0, "c": 1, "d": 0},
        "W3": {"a": 1, "b": 1, "c": 0, "d": 1},
        "W4": {"a": 1, "b": 0, "c": None, "d": 1},
    }
    matrix = tradition.content_characters(coverages)
    assert matrix.taxa == ("W1", "W2", "W3", "W4")
    assert matrix.characters == ("b", "d")
    assert matrix.as_dict() == {"b": [1, 0, 1, 0], "d": [0, 0, 1, 1]}
    assert tradition.content_characters(coverages, min_each=1).characters == (
        "b",
        "c",
        "d",
    )


def test_content_characters_universe_order():
    # Coverages with different items: without `universe`, columns follow the
    # order of first appearance; with it, the given order
    coverages = {
        "W1": {"b": 1, "a": 1},
        "W2": {"a": 0, "b": 0, "c": 1},
        "W3": {"c": 0, "a": 1, "b": 1},
        "W4": {"a": 0, "b": 0, "c": 1},
    }
    assert tradition.content_characters(coverages).characters == ("b", "a")
    matrix = tradition.content_characters(coverages, universe="abc")
    assert matrix.characters == ("a", "b")
    assert matrix.as_dict() == {"a": [1, 0, 1, 0], "b": [1, 0, 1, 0]}

    # Items outside the universe are ignored
    assert tradition.content_characters(coverages, universe="a").characters == ("a",)


def test_content_characters_names_and_missing():
    # Items are named with str(); an item missing from a coverage is missing
    coverages = {
        "W1": {1: 1, (2, 3): 0},
        "W2": {1: 1, (2, 3): 1},
        "W3": {1: 0, (2, 3): 1},
        "W4": {1: 0},
    }
    matrix = tradition.content_characters(coverages, min_each=1)
    assert matrix.characters == ("1", "(2, 3)")
    assert matrix.column("(2, 3)") == [0, 1, 1, None]
    assert tradition.content_characters({}).shape == (0, 0)


def test_adjacencies():
    assert tradition.adjacencies("abXcab", "abc") == {("a", "b"), ("b", "c")}
    assert tradition.adjacencies("a", "a") == set()
    assert tradition.adjacencies("", "abc") == set()
    assert tradition.adjacencies("abc", {"a", "c"}) == {("a", "c")}


def test_adjacency_characters():
    witnesses = {"W1": "abcd", "W2": "abcd", "W3": "acbd", "W4": "acbd", "W5": "ad"}
    matrix = tradition.adjacency_characters(witnesses, "abcd")
    assert matrix.characters == ("a>b", "a>c", "b>c", "b>d", "c>b", "c>d")
    assert matrix.column("a>b") == [1, 1, 0, 0, None]
    # "a>d" is adjacent in one witness only
    assert "a>d" not in matrix.characters
    # With a single witness per state, everything is kept
    matrix = tradition.adjacency_characters(witnesses, "abcd", min_each=1)
    assert matrix.column("a>d") == [0, 0, 0, 0, 1]
    # Items outside the universe are skipped: "d" is dropped
    matrix = tradition.adjacency_characters(witnesses, "abc")
    assert matrix.characters == ("a>b", "a>c", "b>c", "c>b")
    assert matrix.column("a>c") == [0, 0, 1, 1, None]


def test_adjacency_characters_parts():
    # W3 is made of two interleaved parts: "a" and "c" are in different
    # parts, so "a>c" is missing, not absent, and "b>d" is present
    witnesses = {"W1": "abcd", "W2": "abcd", "W3": "abdc", "W4": "acdb"}
    parts_of = {"W3": ["bd", "ac"]}
    matrix = tradition.adjacency_characters(
        witnesses, "abcd", min_each=1, parts_of=parts_of
    )
    assert matrix.column("a>b") == [1, 1, None, 0]
    assert matrix.column("a>c") == [0, 0, 1, 1]
    assert matrix.column("b>c") == [1, 1, None, 0]
    assert matrix.column("b>d") == [0, 0, 1, 0]
    # No witness has both "c" and "d" without "c>d": not informative
    assert "c>d" not in matrix.characters
    matrix = tradition.adjacency_characters(witnesses, "abcd", min_each=1)
    assert matrix.column("a>b") == [1, 1, 1, 0]


def test_adjacency_characters_order_and_names():
    # Columns are sorted by pair, not by name: ("a", "b") < ("a-", "b"),
    # but "a->b" < "a>b"
    witnesses = {
        "W1": ["a", "b", "a-"],
        "W2": ["a", "b", "a-"],
        "W3": ["a-", "b", "a"],
        "W4": ["a-", "b", "a"],
    }
    matrix = tradition.adjacency_characters(witnesses, ["a", "a-", "b"])
    assert matrix.characters == ("a>b", "a->b", "b>a", "b>a-")
    assert matrix.column("a->b") == [0, 0, 1, 1]
    # A key for items that cannot be compared, and another separator
    witnesses = {"W1": [1, "x"], "W2": [1, "x"], "W3": ["x", 1], "W4": ["x", 1]}
    matrix = tradition.adjacency_characters(
        witnesses, [1, "x"], key=str, separator="->"
    )
    assert matrix.characters == ("1->x", "x->1")
    with pytest.raises(TypeError):
        tradition.adjacency_characters(witnesses, [1, "x"])
    assert tradition.adjacency_characters({}, "abc").shape == (0, 0)


def test_restrict():
    matrix = CharacterMatrix(
        ("W1", "W2", "W3", "W4", "W5"),
        ("a", "b"),
        ((1, 1), (1, 0), (0, 1), (0, 0), (0, 1)),
    )
    restricted = tradition.restrict(matrix, ["W4", "W3", "W2", "W1"])
    assert restricted.taxa == ("W4", "W3", "W2", "W1")
    assert restricted.cells == ((0, 0), (0, 1), (1, 0), (1, 1))
    assert tradition.restrict(matrix, ["W1", "W2", "W3", "W5"]).characters == ("a",)
    assert tradition.restrict(matrix, ["W1", "W2"]).shape == (2, 0)
    assert tradition.restrict(matrix, ["W1", "W2"], min_each=1).characters == ("b",)


def test_concat_characters():
    m1 = CharacterMatrix(("W1", "W2"), ("a",), ((1,), (0,)))
    m2 = CharacterMatrix(("W1", "W2"), ("a>b", "b>c"), ((0, 1), (None, 1)))
    combined = tradition.concat_characters([m1, m2])
    assert combined.characters == ("a", "a>b", "b>c")
    assert combined.cells == ((1, 0, 1), (0, None, 1))
    assert tradition.concat_characters([]).shape == (0, 0)
    with pytest.raises(ValueError):
        tradition.concat_characters([m1, CharacterMatrix(("W2", "W1"), (), ((), ()))])


def test_character_blocks():
    position = {item: pos for pos, item in enumerate("abcdef")}
    names = ["a", "c>d", "e>a", "f"]
    assert tradition.character_blocks(names, position, size=2) == [0, 1, 2, 2]
    assert tradition.character_blocks(names, position, size=4) == [0, 0, 1, 1]
    label_of = dict(zip("abcdef", ["I", "I", "I", "J", "J", "J"]))
    assert tradition.character_blocks(names, label_of=label_of) == ["I", "I", "J", "J"]
    # Keys are converted with str(), as in character names
    position = {1: 0, 2: 1, 3: 2}
    assert tradition.character_blocks(["1>2", "3"], position, size=2) == [0, 1]
    # Content characters whose items contain the separator
    position = {"a>b": 0, "c": 5}
    assert tradition.character_blocks(
        ["a>b", "c"], position, size=5, separator=None
    ) == [0, 1]
    with pytest.raises(ValueError):
        tradition.character_blocks(names, position)
    with pytest.raises(ValueError):
        tradition.character_blocks(names, position, size=2, label_of=label_of)
    with pytest.raises(ValueError):
        tradition.character_blocks(names, size=2)


# Export


def test_taxon_name():
    assert tradition.taxon_name("Krka 4 A+E+G") == "Krka_4_A-E-G"
    assert tradition.taxon_name("W1") == "W1"
    assert tradition.taxon_name("a (b)", [(" ", "_"), ("(", ""), (")", "")]) == "a_b"
    assert tradition.taxon_name("a b", []) == "a b"


def test_to_phylip():
    matrix = CharacterMatrix(("W 1", "W+2"), ("a", "b", "c"), ((1, 0, None), (0, 1, 1)))
    assert tradition.to_phylip(matrix) == "2 3\nW_1 10?\nW-2 011\n"
    assert tradition.to_phylip(matrix, rename=str.lower) == "2 3\nw 1 10?\nw+2 011\n"
    assert tradition.to_phylip(CharacterMatrix((), (), ())) == "0 0\n"


def test_to_nexus():
    matrix = CharacterMatrix(("W1", "Witness 2"), ("a", "b"), ((1, None), (0, 1)))
    expected = (
        "#NEXUS\n"
        "begin data;\n"
        "  dimensions ntax=2 nchar=2;\n"
        "  format datatype=restriction missing=? gap=-;\n"
        "  matrix\n"
        "    W1         1?\n"
        "    Witness_2  01\n"
        "  ;\n"
        "end;\n"
    )
    assert tradition.to_nexus(matrix) == expected
    assert tradition.to_nexus(matrix, {}) == expected
    with_sets = tradition.to_nexus(matrix, {"content": (1, 1), "adjacency": (2, 2)})
    assert with_sets == expected + (
        "begin sets;\n  charset content = 1-1;\n  charset adjacency = 2-2;\nend;\n"
    )
    assert tradition.to_nexus(CharacterMatrix((), (), ())).count("\n") == 7


def test_partitions():
    charsets = tradition.charset_ranges({"content": 3, "adjacency": 2})
    assert charsets == {"content": (1, 3), "adjacency": (4, 5)}
    assert tradition.charset_ranges({}) == {}
    assert tradition.partition_nexus(charsets) == (
        "#nexus\nbegin sets;\n  charset content = 1-3;\n  charset adjacency = 4-5;\nend;\n"
    )


def test_mrbayes():
    ages = {"W1": (500, 500), "W2": (300, 400), "W3": (1.5, 2.0)}
    assert tradition.mrbayes_calibrations(ages) == [
        "  calibrate W1 = fixed(500);",
        "  calibrate W2 = uniform(300,400);",
        "  calibrate W3 = uniform(1.5,2.0);",
    ]
    assert tradition.mrbayes_calibrations({}) == []
    groups = {"west": (["W1", "W2"], (600, 900)), "east": (["W3"], None)}
    assert tradition.mrbayes_constraints(groups) == [
        "  constraint west = W1 W2;",
        "  constraint east = W3;",
        "  calibrate west = uniform(600,900);",
        "  prset topologypr=constraints(west,east);",
    ]
    assert tradition.mrbayes_constraints({}) == []


# Determinism

SCRIPT = """
import random
from seqsim import tradition

rng = random.Random(1)
pool = [f"s{i:02d}" for i in range(40)]
witnesses = {}
for idx in range(10):
    seq = [s for s in pool if rng.random() > 0.2]
    i, j = sorted(rng.sample(range(len(seq)), 2))
    seq = seq[:i] + seq[j:] + seq[i:j]
    witnesses[f"W {idx}"] = seq
frame = tradition.consensus_order(witnesses.values())
label_of = {item: f"C{int(item[1:]) // 10}" for item in frame}
coverages = {name: tradition.coverage(seq, frame, label_of) for name, seq in witnesses.items()}
parts = {"W 0": [witnesses["W 0"][:15], witnesses["W 0"][15:]]}
content = tradition.content_characters(coverages)
adjacency = tradition.adjacency_characters(witnesses, frame, parts_of=parts)
print(frame)
print(sorted(tradition.covered_labels(witnesses["W 1"], frame, label_of)))
print(tradition.merge_parts([("A", pool[:10]), ("B", pool[10:30]), ("C", pool[5:12])]))
print(tradition.to_phylip(tradition.concat_characters([content, adjacency])))
print(tradition.character_blocks(adjacency.characters, label_of=label_of))
"""


def test_independent_of_hash_seed():
    root = pathlib.Path(__file__).parent.parent
    outputs = []
    for seed in ("0", "1", "12345"):
        env = dict(os.environ, PYTHONHASHSEED=seed)
        env["PYTHONPATH"] = os.pathsep.join(
            [str(root / "src"), env.get("PYTHONPATH", "")]
        )
        result = subprocess.run(
            [sys.executable, "-c", SCRIPT],
            env=env,
            capture_output=True,
            text=True,
            check=True,
        )
        outputs.append(result.stdout)
    assert outputs[0] == outputs[1] == outputs[2]
    assert "W_9 " in outputs[0]
