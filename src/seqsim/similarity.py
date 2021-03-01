"""
Module implementing different similarity scores.
"""

# Import Python standard libraries
import difflib
from typing import Hashable, Sequence

# Import from local modules
from .common import _nwise, _indices

# TODO: allow normalized version of similarities, in range [0..1]


def fast_birnbaum_simil(seq_x: Sequence[Hashable], seq_y: Sequence[Hashable]) -> int:
    """
    Compute the Birnbaum similarity score with the fast method.

    This implementation uses the experimental method we developed following
    the description in Birnbaum (2003), which is much faster and less
    memory-intensive than the one implemented in the `birnbaum_simil()`
    function. While in most cases the results are comparable, and the ones
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

    # Make sure `seq_x` is shorter or equal in length to `seq_y`
    if len(seq_x) < len(seq_y):
        seq_x, seq_y = seq_y, seq_x

    # Get the opcodes from `SequenceMatcher` and drop the last one (a dummy)
    sm = difflib.SequenceMatcher(None, seq_x, seq_y)
    blocks = sm.get_matching_blocks()[:-1]

    # Sum sizes and return
    sizes = [match.size for match in blocks]
    formula = [(v * (v + 1)) // 2 for v in sizes]

    return sum(formula)


def birnbaum_simil(seq_x: Sequence[Hashable], seq_y: Sequence[Hashable]) -> int:
    """
    Compute the Birnbaum similarity score with the original method.

    This implementation uses the method first developed for Göransson et al.,
    following the description in Birnbaum (2003). Note that, on large
    sequences, it is much slower and more memory-intensive than the one
    implemented in the `fast_birnbaum_simil` method.  While in most cases
    the results are comparable, the values are *not* identical.

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

    # Make sure `seq_x` is shorter or equal in length to `seq_y`
    if len(seq_x) < len(seq_y):
        seq_x, seq_y = seq_y, seq_x

    # Cast to lists to make sure it works with arbitrary hashable elements
    seq_x = [element for element in seq_x]
    seq_y = [element for element in seq_y]

    similarity: int = 0
    subseq_len = len(seq_y)
    while subseq_len:
        # iterate though all subsequences of N
        for Nsubseq in _nwise(seq_y, subseq_len):
            # test whether each of these subsequences matches a section of M
            # (go to each instance of Nsubseq[0] in M and test the
            # appropriate length slice from there)
            for i in _indices(seq_x, Nsubseq[0]):
                if seq_x[i : i + subseq_len] == Nsubseq:
                    similarity += 1
        subseq_len -= 1

    return similarity
