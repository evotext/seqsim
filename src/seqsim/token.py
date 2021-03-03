"""
Module implementing various methods for similarity and distance from token methods.

Most of the methods are commonly used in string comparison, such as Jaccard
index, but in this module we need to make sure we can operate on arbitrary
iterable data structures.
"""

# Import Python standard libraries
from typing import Hashable, Sequence
import logging

# Import 3rd-party libraries
import textdistance

# Import local modules
from .common import collect_subseqs

# TODO: have a jaccard similarity?
def jaccard_dist(
    seq_x: Sequence[Hashable], seq_y: Sequence[Hashable], normal: bool = False
) -> float:
    """
    Computes a Jaccard distance between two sequences.

    The function accepts the `normal` parameter to have calls equivalent to those
    of other methods, but it is redundant as the Jaccard distance is already
    in range [0..1].

    :param seq_x: The first sequence to be compared.
    :param seq_y: The second sequence to be compared.
    :param normal: Dummy parameter, see comment above.
    :return: The Jaccard distance between the two sequences.
    """

    intersection = len(set(seq_x).intersection(seq_y))
    union = len(seq_x) + len(seq_y) - intersection

    if normal:
        logging.warning(
            "Jaccard distance is always in [0..1] range, no need for `normal` parameter."
        )

    return 1.0 - (float(intersection) / union)


# TODO: rename appropriately with other methods, consider ngram usage
# TODO: have a subseq_jaccard similarity?
def subseq_jaccard_dist(
    seq_x: Sequence[Hashable], seq_y: Sequence[Hashable], normal: bool = False
) -> float:
    """
    Computes a Jaccard distance between two sequences using sub-sequence occurrence.

    The function accepts the `normal` parameter to have calls equivalent to those
    of other methods, but it is redundant as the Jaccard distance is already
    in range [0..1].

    :param seq_x: The first sequence to be compared.
    :param seq_y: The second sequence to be compared.
    :param normal: Dummy parameter, see comment above.
    :return: The Subseq-Jaccard distance between the two sequences.
    """

    subseqs1 = collect_subseqs(seq_x)
    subseqs2 = collect_subseqs(seq_y)

    # From the longest subseq, which is the length of the longest sequence,
    # collect all subsequences of that given length in both sets, compute the
    # Jaccard index, correct it by length of the subsequence (so that longer ones
    # will score higher) and update the internal results before returning.
    # TODO: collect beforehand in a list, so we don´t repeat the comprehension
    # TODO: convert to tuple beforehand as well?
    jaccard_scores = []
    max_length = max([len(seq_x), len(seq_y)])
    for length in range(max_length, 0, -1):
        l_subseq1 = [tuple(ss) for ss in subseqs1 if len(ss) == length]
        l_subseq2 = [tuple(ss) for ss in subseqs2 if len(ss) == length]

        intersection = len(set(l_subseq1).intersection(l_subseq2))
        union = len(l_subseq1) + len(l_subseq2) - intersection
        jaccard_scores.append((float(intersection) / union) * length)

    # Compute the denominator, as the highest possible value
    den = (max_length * (max_length + 1)) / 2.0

    if normal:
        logging.warning(
            "Subseq-Jaccard distance is always in [0..1] range, no need for `normal` parameter."
        )

    return (1.0 - (sum(jaccard_scores) / den)) ** max_length


def sorensen_dist(
    seq_x: Sequence[Hashable], seq_y: Sequence[Hashable], normal: bool = False
) -> float:
    """
    Computes a distance between two sequences based on the Sørensen–Dice coefficient.

    The function accepts the `normal` parameter to have calls equivalent to those
    of other methods, but it is redundant as the Sørensen–Dice distance is already
    in range [0..1].

    See: https://en.wikipedia.org/wiki/S%C3%B8rensen%E2%80%93Dice_coefficient

    :param seq_x: The first sequence to be compared.
    :param seq_y: The second sequence to be compared.
    :param normal: Dummy parameter, see comment above.
    :return: The Sørensen–Dice distance between the two sequences.
    """

    if normal:
        logging.warning(
            "Sørensen–Dice distance is always in [0..1] range, no need for `normal` parameter."
        )

    return 1.0 - textdistance.Sorensen(external=False)(seq_x, seq_y)
