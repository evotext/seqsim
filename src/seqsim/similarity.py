"""
Module implementing different similarity scores.
"""

# Import Python standard libraries
import difflib
from typing import Hashable, Sequence

# TODO: allow normalized version of similarities, in range [0..1]
# TODO: test symmetries


def fast_birnbaum_simil(seq_x: Sequence[Hashable], seq_y: Sequence[Hashable]) -> int:
    """
    Compute the Birnbaum similarity score with the fast method.

    This implementation uses the experimental method we developed following
    the description in Birnbaum (2003), which is much faster and less
    memory-intensive than the one implemented in the `birnbaum_simil()`
    method. While in most cases the results are comparable, and the ones
    provided by this method might be considered more adequate due to their
    handling of duplicate information, the values are *not* identical.

    See: Birnbaum, David J. (2003). "Computer-Assisted Analysis and
        Study of the Structure of Mixed-Content Miscellanies".
        Scripta & Scripta 1:15-64.

    @param seq_x: The first sequence to be compared.
    @param seq_y: The second sequence to be compared.
    @return: The similarity score between the two sequences. The higher
        the similarity score, the more similar the two sequences are;
        a similarity score of zero is the theoretical maximum difference
        between two sequences.
    """

    # If the sequences are equal, we can just compute the score from length
    if seq_x == seq_y:
        length = len(seq_x)
        return (length * (length + 1)) // 2

    # Get the opcodes from `SequenceMatcher` and drop the last one (a dummy)
    sm = difflib.SequenceMatcher(None, seq_x, seq_y)
    blocks = sm.get_matching_blocks()[:-1]

    # Sum sizes and return
    sizes = [match.size for match in blocks]
    formula = [(v * (v + 1)) // 2 for v in sizes]

    return sum(formula)
