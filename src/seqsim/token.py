"""
Module implementing methods for sequence dissimilarity based on tokens.

These methods compare the elements (or sub-sequences) that two sequences
share, such as the Jaccard index, operating on arbitrary sequences of hashable
elements. See the `edit` module for the naming convention of the functions.
"""

# Import Python standard libraries
from collections import Counter
from typing import Hashable, Sequence

# Import local modules
from .common import empty_dissim, equivalent_string


def jaccard_dissim(
    seq_x: Sequence[Hashable], seq_y: Sequence[Hashable], *, normal: bool = False
) -> float:
    """
    Computes the Jaccard dissimilarity between two sequences.

    The dissimilarity is one minus the Jaccard index of the sets of elements
    of both sequences. While the Jaccard distance is a true distance on sets,
    order and repetition are ignored, so different sequences can have a
    dissimilarity of zero (e.g., `"ab"` and `"ba"`, or `"a"` and `"aa"`).

    Results are always in range [0..1], so `normal` has no effect.

    Example
    ********

    .. code-block:: python

        >>> seqsim.token.jaccard_dissim("abc", "bcde")
        0.6

    References
    ***********

    Tan PN, Steinbach M, Kumar V (2005). Introduction to Data Mining. ISBN 0-321-32136-7.

    :param seq_x: The first sequence to be compared.
    :param seq_y: The second sequence to be compared.
    :param normal: Ignored, as results are always in range [0..1].
    :return: The Jaccard dissimilarity between the two sequences.
    """

    empty = empty_dissim(seq_x, seq_y)
    if empty is not None:
        return empty

    set_x, set_y = set(seq_x), set(seq_y)

    return 1.0 - (len(set_x & set_y) / len(set_x | set_y))


def subseq_jaccard_dissim(
    seq_x: Sequence[Hashable], seq_y: Sequence[Hashable], *, normal: bool = False
) -> float:
    """
    Computes a Jaccard dissimilarity between two sequences using sub-sequences.

    For each length `n` from 1 to the length of the longest sequence, the
    Jaccard index is computed on the multisets of contiguous sub-sequences of
    length `n` of both sequences. The similarity is the mean of these indices,
    weighted by `n` so that longer shared sub-sequences count more, and the
    dissimilarity is one minus this similarity. Identical sequences, and only
    identical sequences, have a dissimilarity of zero; the measure does not
    satisfy the triangle inequality.

    Results are always in range [0..1], so `normal` has no effect.

    Example
    ********

    .. code-block:: python

        >>> seqsim.token.subseq_jaccard_dissim("abc", "bcde")
        0.91

    References
    ***********

    Tan PN, Steinbach M, Kumar V (2005). Introduction to Data Mining. ISBN 0-321-32136-7.

    :param seq_x: The first sequence to be compared.
    :param seq_y: The second sequence to be compared.
    :param normal: Ignored, as results are always in range [0..1].
    :return: The Subseq-Jaccard dissimilarity between the two sequences.
    """

    empty = empty_dissim(seq_x, seq_y)
    if empty is not None:
        return empty

    # Strings are much faster to slice and hash than tuples
    str_x, str_y = equivalent_string(seq_x, seq_y)
    max_length = max(len(str_x), len(str_y))

    weighted_sum = 0.0
    for length in range(1, max_length + 1):
        counter_x = Counter(
            str_x[i : i + length] for i in range(len(str_x) - length + 1)
        )
        counter_y = Counter(
            str_y[i : i + length] for i in range(len(str_y) - length + 1)
        )

        # Use multisets for both the intersection and the union, so that
        # repeated sub-sequences are counted consistently
        intersection = sum((counter_x & counter_y).values())
        union = sum((counter_x | counter_y).values())
        weighted_sum += length * (intersection / union)

    # The highest possible value is the sum of all weights
    return 1.0 - (weighted_sum / ((max_length * (max_length + 1)) / 2))


def sorensen_dissim(
    seq_x: Sequence[Hashable], seq_y: Sequence[Hashable], *, normal: bool = False
) -> float:
    """
    Computes a dissimilarity between two sequences based on the Sørensen–Dice coefficient.

    The dissimilarity is one minus the Sørensen–Dice coefficient of the
    multisets of elements of both sequences. Order is ignored, so different
    sequences can have a dissimilarity of zero (e.g., `"ab"` and `"ba"`), and
    the measure does not satisfy the triangle inequality.

    Results are always in range [0..1], so `normal` has no effect.

    Example
    ********

    .. code-block:: python

        >>> seqsim.token.sorensen_dissim("abc", "bcde")
        0.4285714285714286

    References
    ***********

    Kondrak, Grzegorz; Marcu, Daniel; Knight, Kevin (2003). "Cognates Can Improve Statistical
    Translation Models" (PDF). Proceedings of HLT-NAACL 2003: Human Language Technology
    Conference of the North American Chapter of the Association for Computational
    Linguistics. pp. 46–48.

    Sørensen, T. (1948). "A method of establishing groups of equal amplitude in plant
    sociology based on similarity of species and its application to analyses of the
    vegetation on Danish commons". Kongelige Danske Videnskabernes Selskab. 5 (4): 1–34.

    :param seq_x: The first sequence to be compared.
    :param seq_y: The second sequence to be compared.
    :param normal: Ignored, as results are always in range [0..1].
    :return: The Sørensen–Dice dissimilarity between the two sequences.
    """

    empty = empty_dissim(seq_x, seq_y)
    if empty is not None:
        return empty

    intersection = sum((Counter(seq_x) & Counter(seq_y)).values())

    return 1.0 - (2.0 * intersection / (len(seq_x) + len(seq_y)))
