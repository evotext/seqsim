"""
Module implementing various methods for similarity and distance from sequence methods.

Most of the methods are commonly used in string comparison, such as
RatcliffObershelp, but in this module we need to make sure we can operate on arbitrary
iterable data structures.
"""

# Import Python standard libraries
from typing import Hashable, Sequence
import logging

# Import 3rd-party libraries
import textdistance

# TODO: add a pure ratcliff_obershelp similarity?
# TODO: multiple sequences?
def ratcliff_obershelp(
    seq_x: Sequence[Hashable], seq_y: Sequence[Hashable], normal: bool = False
) -> float:
    """
    Computes a distance between two sequences based on the Ratcliff-Obershelp similarity.

    The function accepts the `normal` parameter to have calls equivalent to those
    of other methods, but it is redundant as the Ratcliff-Obershelp distance is already
    in range [0..1].

    See: https://en.wikipedia.org/wiki/Gestalt_Pattern_Matching

    :param seq_x: The first sequence to be compared.
    :param seq_y: The second sequence to be compared.
    :param normal: Dummy parameter, see comment above.
    :return: The Ratcliff-Obershelp distance between the two sequences.
    """

    # As the method uses .find, we need strings
    # TODO: rewrite proper solution generalizing for all sequences
    seq_x = "".join([str(elem) for elem in seq_x])
    seq_y = "".join([str(elem) for elem in seq_y])

    if normal:
        logging.warning(
            "Sørensen–Dice distance is always in [0..1] range, no need for `normal` parameter."
        )

    return 1.0 - textdistance.RatcliffObershelp(external=False)(seq_x, seq_y)
