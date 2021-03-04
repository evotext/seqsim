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
    "seq_x,seq_y,expected_x,expected_y",
    [
        ["kitten", "sitting", "kitten", "sitting"],
        ["kitten", [c for c in "sitting"], "326604", "5266241"],
        ["kitten", ["si", "tt", "ing"], "316604", "572"],
        [(1, 2, 3), [1, 2, 3], "012", "012"],
        [(1, 2, 3, 4, 5), (1, 2, 4, 3, 6, 7), "01234", "013256"],
        [(1, 2, 3), ["a", "b", "c", "d"], "012", "3456"],
    ],
)
def test_equivalent_string(seq_x, seq_y, expected_x, expected_y):
    eq_x, eq_y = seqsim.common.equivalent_string(seq_x, seq_y)
    assert eq_x == expected_x
    assert eq_y == expected_y


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
