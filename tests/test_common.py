"""
test_common
===========

Tests for the `common` module of the `seqsim` package.
"""

# TODO: add empty strings

# Import Python standard libraries
import random
import pytest

# Import the library being tested
import seqsim


@pytest.mark.parametrize(
    "seq_x,seq_y",
    [
        ["kitten", [c for c in "sitting"]],
        ["kitten", ["si", "tt", "ing"]],
        [(1, 2, 3), [1, 2, 3]],
        [(1, 2, 3, 4, 5), (1, 2, 4, 3, 6, 7)],
        [(1, 2, 3), ["a", "b", "c", "d"]],
        [[1, 23], [12, 3]],
        [["ab"], ["a", "b"]],
        [[None], ["None"]],
        [["\t", "\n", " "], ["\t", " "]],
    ],
)
def test_equivalent_string(seq_x, seq_y):
    eq_x, eq_y = seqsim.common.equivalent_string(seq_x, seq_y)

    # Lengths are preserved
    assert len(eq_x) == len(seq_x)
    assert len(eq_y) == len(seq_y)

    # Two positions share a character if and only if the elements are equal
    elems = list(seq_x) + list(seq_y)
    chars = eq_x + eq_y
    for i, elem_i in enumerate(elems):
        for j, elem_j in enumerate(elems):
            assert (elem_i == elem_j) == (chars[i] == chars[j])

    # No printable or whitespace characters are used for the mapping
    assert not any(c.isprintable() or c.isspace() for c in chars)


def test_equivalent_string_strings():
    assert seqsim.common.equivalent_string("kitten", "sitting") == (
        "kitten",
        "sitting",
    )


def test_equivalent_string_long():
    """
    Test `equivalent_string()` as above, but with a big set of elements.
    """

    seq_x = list(range(0, 1000))
    seq_y = list(range(500, 2500))

    random.seed("seqsim")
    random.shuffle(seq_x)
    random.shuffle(seq_y)

    eq_x, eq_y = seqsim.common.equivalent_string(seq_x, seq_y)

    assert len(eq_x) == 1000
    assert len(eq_y) == 2000
    assert len(set(eq_x + eq_y)) == 2500
