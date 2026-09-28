"""
test_compression
================

Tests for the `compression` module of the `seqsim` package.
"""

# Import Python standard libraries
import random

import pytest

# Import the library being tested
from seqsim import compression


@pytest.mark.parametrize(
    "seq_x,seq_y",
    [
        ["kitten", "sitting"],
        [(1, 2, 3, 4, 5), (1, 2, 4, 3, 6, 7)],
        [(1, 2, 3), ["a", "b", "c", "d"]],
    ],
)
def test_lzma_ncd_symmetry(seq_x, seq_y):
    assert compression.lzma_ncd_dissim(seq_x, seq_y) == compression.lzma_ncd_dissim(
        seq_y, seq_x
    )


def test_lzma_ncd_long_sequences():
    # For sequences long enough for the compressor to exploit repetition,
    # identical sequences are close to zero and more different sequences
    # score higher
    rng = random.Random("seqsim")
    base = [rng.randrange(8) for _ in range(300)]
    other = [rng.randrange(8) for _ in range(300)]
    half = base[:150] + other[150:]

    same = compression.lzma_ncd_dissim(base, list(base))
    partial = compression.lzma_ncd_dissim(base, half)
    different = compression.lzma_ncd_dissim(base, other)

    assert same < 0.1
    assert same < partial < different


def test_lzma_ncd_arbitrary_elements():
    # Sequences of elements that are not characters, including more than
    # 256 distinct elements
    seq_x = list(range(1000))
    seq_y = list(range(500, 1500))
    assert 0.0 < compression.lzma_ncd_dissim(seq_x, seq_y) <= 1.1
    assert compression.lzma_ncd_dissim([None, (1, 2)], [None, (1, 2)]) >= 0.0


def test_lzma_ncd_normal_clipped():
    # Two very short sequences can have a raw NCD above 1.0
    assert 0.0 <= compression.lzma_ncd_dissim("a", "b", normal=True) <= 1.0


@pytest.mark.parametrize(
    "seq_x,seq_y,expected,tol",
    [
        ["kitten", "sitting", 0.101340, 1e-6],
        [(1, 2, 3), [1, 2, 3], 0.0, 0.0],
        [(1, 2, 3, 4, 5), (1, 2, 4, 3, 6, 7), 0.114430, 1e-6],
        [(1, 2, 3), ["a", "b", "c", "d"], 0.407464, 1e-6],
    ],
)
def test_entropy_ncd(seq_x, seq_y, expected, tol):
    # Values match `textdistance` 4.5, from which the method was ported
    assert compression.entropy_ncd_dissim(seq_x, seq_y) == pytest.approx(
        expected, abs=tol
    )
