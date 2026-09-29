"""
Module for defining common functions and variables used in different circumstances.

This module holds functions that are used by different methods, such as the
mapping of sequences of arbitrary hashable elements to equivalent strings,
along with more low-level utilities for handling sub-sequences.
"""

# Import Python standard libraries
from typing import Hashable, Iterator, Optional, Sequence, Tuple
import itertools


def empty_dissim(
    seq_x: Sequence[Hashable], seq_y: Sequence[Hashable]
) -> Optional[float]:
    """
    Returns the dissimilarity for comparisons involving empty sequences.

    All dissimilarities in range [0..1] follow the same convention: two empty
    sequences are identical (0.0), and an empty sequence is maximally
    different from any non-empty sequence (1.0). If neither sequence is
    empty, `None` is returned.

    :param seq_x: The first sequence to be compared.
    :param seq_y: The second sequence to be compared.
    :return: The dissimilarity, or `None` if neither sequence is empty.
    """

    empty_x, empty_y = len(seq_x) == 0, len(seq_y) == 0
    if empty_x and empty_y:
        return 0.0
    if empty_x or empty_y:
        return 1.0

    return None


# Unicode Private Use Areas, used for mapping arbitrary hashable elements to
# single characters that cannot collide with any "real" (assigned) character
_PRIVATE_USE_RANGES = ((0xE000, 0xF8FF), (0xF0000, 0xFFFFD), (0x100000, 0x10FFFD))
_MAX_EQUIVALENT_SYMBOLS = sum(end - start + 1 for start, end in _PRIVATE_USE_RANGES)


def _private_use_chars() -> Iterator[str]:
    """
    Yields all the characters in the Unicode Private Use Areas, in order.
    """

    for start, end in _PRIVATE_USE_RANGES:
        for codepoint in range(start, end + 1):
            yield chr(codepoint)


def equivalent_string(
    seq_x: Sequence[Hashable], seq_y: Sequence[Hashable]
) -> Tuple[str, str]:
    """
    Returns a string equivalent to a sequence, for comparison.

    Some methods are most efficiently implemented on strings, while `seqsim` is
    designed to offer all methods of comparison for generic sequences of hashable
    elements, so in some cases it is necessary to convert a sequence to an
    equivalent string. Using a normal `str` conversion is not possible or
    satisfactory for a number of reasons, including elements not having a string
    representation, and individual string representations of different lengths
    and potentially overlapping (consider cases like `[1, 12, 123, 23]`).

    This function accepts a pair of sequences and returns an equivalent
    textual representation, that is, a pair of strings where the order is
    preserved and each element is mapped to a single, unique character. Two
    elements are mapped to the same character if and only if they are equal
    (following Python's `==` and `hash()` semantics). Characters are taken from
    the Unicode Private Use Areas, assigned in the order of the elements sorted
    by type name and `repr()`, so the mapping is deterministic and does not
    depend on the order of the arguments.

    If two strings are passed, the same strings will be returned. Note
    that in case of mixed types (such as a string and a list of
    strings), strings will be considered sequences of characters
    and will be modified upon return, as the tokens of the
    second sequence could be of length over one character (e.g.,
    `"abc"` and `["a", "bc"]`).

    Example
    ********

    .. code-block:: python

        >>> x, y = seqsim.common.equivalent_string([1, 2, 3], [1, 2, 4, 5])
        >>> [ord(c) for c in x], [ord(c) for c in y]
        ([57344, 57345, 57346], [57344, 57345, 57347, 57348])

    :param seq_x: The first sequence to be mapped to an equivalent
        string.
    :param seq_y: The second sequence to be mapped to an equivalent
        string.
    :return: A tuple of two strings equivalent, for matters of
        comparison and distance computation, to the provided
        sequences.
    :raises ValueError: If the sequences have more distinct elements than
        there are characters in the Unicode Private Use Areas.
    """

    # Don't need to apply to strings
    if isinstance(seq_x, str) and isinstance(seq_y, str):
        return seq_x, seq_y

    # Collect the distinct elements, sorted by type and representation so that
    # the mapping does not depend on the order of the arguments (ties, which
    # are very unlikely, are kept in order of first appearance)
    elements = sorted(
        dict.fromkeys(itertools.chain(seq_x, seq_y)),
        key=lambda element: (type(element).__qualname__, repr(element)),
    )
    if len(elements) > _MAX_EQUIVALENT_SYMBOLS:
        raise ValueError(
            f"Cannot map {len(elements)} distinct elements to single characters "
            f"(maximum is {_MAX_EQUIVALENT_SYMBOLS})."
        )

    mapper = dict(zip(elements, _private_use_chars()))

    return (
        "".join(mapper[element] for element in seq_x),
        "".join(mapper[element] for element in seq_y),
    )


def lcs_length(seq_x: Sequence[Hashable], seq_y: Sequence[Hashable]) -> int:
    """
    Returns the length of the longest common subsequence of two sequences.

    The subsequence does not need to be contiguous. It is computed with the
    standard dynamic programming algorithm in O(len(x) * len(y)) time and
    O(len(y)) memory.

    Example
    ********

    .. code-block:: python

        >>> seqsim.common.lcs_length("kitten", "sitting")
        4

    :param seq_x: The first sequence.
    :param seq_y: The second sequence.
    :return: The length of the longest common subsequence.
    """

    len_y = len(seq_y)
    prev = [0] * (len_y + 1)
    for elem_x in seq_x:
        curr = [0] * (len_y + 1)
        for j, elem_y in enumerate(seq_y, start=1):
            if elem_x == elem_y:
                curr[j] = prev[j - 1] + 1
            else:
                curr[j] = max(prev[j], curr[j - 1])
        prev = curr

    return prev[len_y]
