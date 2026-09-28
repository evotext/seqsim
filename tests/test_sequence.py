"""
test_sequence
=============

Tests for the `sequence` module of the `seqsim` package.
"""

# Import Python standard libraries
import pytest

# Import the library being tested
from seqsim import sequence


@pytest.mark.parametrize(
    "seq_x,seq_y,expected",
    [
        ("kitten", "sitting", 1 - 8 / 13),  # "itt" + "n"
        ((1, 2, 3), [1, 2, 3], 0.0),
        ((1, 2, 3, 4, 5), (1, 2, 4, 3, 6, 7), 1 - 6 / 11),
        ((1, 2, 3), ["a", "b", "c", "d"], 1.0),
        ("abc", "bcde", 1 - 4 / 7),
    ],
)
def test_ratcliff_obershelp(seq_x, seq_y, expected):
    assert sequence.ratcliff_obershelp_dissim(seq_x, seq_y) == pytest.approx(expected)


@pytest.mark.parametrize(
    "seq_x,seq_y",
    [
        [[1, 23], [12, 3]],
        [["ab"], ["a", "b"]],
        [[None], ["None"]],
    ],
)
def test_ratcliff_obershelp_no_str_collisions(seq_x, seq_y):
    # In 0.3.1 elements were joined as strings, making these pairs identical
    assert sequence.ratcliff_obershelp_dissim(seq_x, seq_y) == 1.0


def test_ratcliff_obershelp_long_sequences():
    # The difflib "autojunk" heuristic must not affect long sequences
    seq_x = ["a"] * 150 + ["b"] * 150
    assert sequence.ratcliff_obershelp_dissim(seq_x, list(seq_x)) == 0.0
    assert sequence.ratcliff_obershelp_dissim(seq_x, seq_x[:-1]) < 0.01
