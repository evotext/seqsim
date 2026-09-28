"""
Module for defining common functions and variables used in different circumstances.

This module works as a big repository of all the functions and variables that
are used by different methods (such as for the computation of an edit
distance using the Wagner-Fischer algorithm), including more low-level and
book-keeping functions such as for interfacing with the system.
"""

# Import Python standard libraries
from typing import Hashable, Iterator, List, Optional, Sequence, Tuple
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


# TODO: replace with the ngram collector module
def collect_subseqs(sequence: Sequence, sort: bool = True) -> List[Sequence]:
    """
    Collects all possible sub-sequences in a given sequence.

    When sorting is requested, sub-sequences will first be sorted by their length and,
    later, by comparing one with the other. Mixing types, like strings and integers, can
    lead to unexpected results and is not suggested if the type cannot be guaranteed.

    Note that this function performs simple comprehensions, neither using padding
    symbols nor the more complex methods n-gram collection methods ultimately based on
    `ngram_iter()`.

    Example
    ********

    .. code-block:: python

        >>> seqsim.common.collect_subseqs('abcde')
        ['a', 'b', 'c', 'd', 'e', 'ab', 'bc', 'cd', 'de', 'abc', 'bcd', 'cde', 'abcd', 'bcde', 'abcde']

    :param sequence: The sequence that shall be converted into it's ngram-representation.
    :param sort: Whether to sort the list of ngrams by length and by identity
        (default: True).
    :return: A list of all ngrams of the input sequence.
    """

    # Cache the length of the sequence
    length = len(sequence)

    # Set the starting index
    idx = 0

    # define the output list
    ret = []

    # start the while loop
    while idx != length and idx < length:
        # copy the sequence
        new_sequence = sequence[idx:length]

        # append the sequence to the output list
        ret += [new_sequence]

        # loop over the new sequence
        for j in range(1, len(new_sequence)):
            ret += [new_sequence[:j]]
            ret += [new_sequence[j:]]

        # increment idx and decrement length
        idx += 1
        length -= 1

    if sort:
        # We try to sort normally; if there is a TypeError, such as when the list has mixed
        # ints and strings, we sort by the string representation of all elements
        # TODO: do it in a better way
        try:
            ret = sorted(ret, key=lambda e: (len(e), e))
        except TypeError:
            ret = sorted(ret, key=lambda e: (len(str(e)), str(e)))

    return ret


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


# TODO: properly rewrite, perhaps using equivalent_string()
def sequence_find(hay: Sequence, needle: Sequence) -> Optional[int]:
    """
    Return the index for starting index of a sub-sequence within a sequence.

    The function is intended to work similarly to the built-in `.find()` method for
    Python strings, but accepting all types of sequences (including different types
    for `hay` and `needle`).

    Example
    ********

    .. code-block:: python

        >>> seqsim.common.sequence_find([1, 2, 3, 4, 5], [2, 3])
        1

    :param hay: The sequence to be searched within.
    :param needle: The sub-sequence to be located in the sequence.
    :return: The starting index of the sub-sequence in the sequence, or `None` if not
             found.
    """
    # Cache `needle` length and have it as a tuple already
    len_needle = len(needle)
    t_needle = tuple(needle)

    # Iterate over all sub-lists (or sub-tuples) of the correct length and check
    # for matches
    for i in range(len(hay) - len_needle + 1):
        if tuple(hay[i : i + len_needle]) == t_needle:
            return i

    return None
