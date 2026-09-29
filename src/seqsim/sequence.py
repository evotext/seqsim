"""
Module implementing methods for sequence dissimilarity based on sequence matching.

These methods, such as Ratcliff-Obershelp, operate on arbitrary sequences of
hashable elements. See the `edit` module for the naming convention of the
functions.
"""

# Import Python standard libraries
from typing import Hashable, Sequence
import difflib

# Import local modules
from ._measure import measure


@measure(
    key="ratcliff_obershelp",
    kind="dissim",
    triangle="no",
    triangle_example=("bcc", "baca", "aa"),
    empty="max",
    symmetrize="min",
)
def ratcliff_obershelp_dissim(
    seq_x: Sequence[Hashable], seq_y: Sequence[Hashable]
) -> float:
    """
    Computes a dissimilarity between two sequences based on the Ratcliff-Obershelp similarity.

    The similarity is twice the number of matching elements divided by the
    total number of elements, where matches are found by recursively taking
    the longest common sub-sequence (as in Python's `difflib`). As the
    matching depends on the order of the arguments, it is computed in both
    orders and the highest similarity is used, so that the measure is
    symmetric. It does not satisfy the triangle inequality.

    Example
    ********

    .. code-block:: python

        >>> seqsim.sequence.ratcliff_obershelp_dissim("abc", "bcde")
        0.4285714285714286

    References
    ***********

    John W. Ratcliff and David Metzener: Pattern Matching: The Gestalt Approach, Dr.
    Dobb's Journal, Issue 46, July 1988

    :param seq_x: The first sequence to be compared.
    :param seq_y: The second sequence to be compared.
    :return: The Ratcliff-Obershelp dissimilarity between the two sequences.
    """

    # `SequenceMatcher` operates directly on sequences of hashable elements;
    # the automatic junk heuristic is disabled, as it would otherwise change
    # the results for sequences with 200 or more elements.
    return 1.0 - difflib.SequenceMatcher(None, seq_x, seq_y, autojunk=False).ratio()
