"""
Module implementing various methods for similarity and distance from sequence methods.

Most of the methods are commonly used in string comparison, such as
RatcliffObershelp, but in this module we need to make sure we can operate on arbitrary
iterable data structures.
"""

# Import Python standard libraries
from typing import Hashable, Sequence
import difflib
import logging

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

    Example
    ********

    .. code-block:: python

        >>> seqsim.sequence.ratcliff_obershelp("abc", "bcde")
        0.4285714285714286

    References
    ***********

    John W. Ratcliff and David Metzener: Pattern Matching: The Gestalt Approach, Dr.
    Dobb's Journal, Issue 46, July 1988

    :param seq_x: The first sequence to be compared.
    :param seq_y: The second sequence to be compared.
    :param normal: Dummy parameter, see comment above.
    :return: The Ratcliff-Obershelp distance between the two sequences.
    """

    # `SequenceMatcher` operates directly on sequences of hashable elements;
    # the automatic junk heuristic is disabled, as it would otherwise change
    # the results for sequences with 200 or more elements.
    ratio = difflib.SequenceMatcher(None, seq_x, seq_y, autojunk=False).ratio()

    if normal:
        logging.warning(
            "Ratcliff-Obershelp distance is always in [0..1] range, no need for `normal` parameter."
        )

    return 1.0 - ratio
