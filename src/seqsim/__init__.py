"""
Main module of the `seqsim` package.

We follow the mathematical definitions for distinguishing between "similarity"
and "distance", as the latter must have the following properties:

  * positivity: d(x,y) >= 0
  * symmetry: d(x,y) = d(y,x)
  * identity-discerning: d(x,y) = 0 => x = y
  * triangle inequality: d(x,z) <= d(x,y) + d(y,z)
"""

# Version of the `seqsim` package
__author__ = "Tiago Tresoldi, Luke Maurits, Michael Dunn"
__email__ = "tiago.tresoldi@lingfil.uu.se"
__version__ = "0.3"

# Import Python standard libraries
from typing import Hashable, Sequence
import itertools

# Import local modules
from . import edit
from . import token
from . import sequence
from . import compression
from .ngrams import ngrams_iter, get_all_ngrams_by_order

# List all methods available and the functions that offers them;
# the dictionary is used by the `distance` function below, but it
# is also a convenient single point of reference for users. When
# a method is available both a similarity and as a distance measure,
# we favor the latter in the label (shorter).
METHODS = {
    "arith_ncd": compression.arith_ncd,
    "birnbaum_simil": edit.birnbaum_simil,
    "birnbaum": edit.birnbaum_dist,
    "bulk_delete": edit.bulk_delete_dist,
    "damerau": edit.levdamerau_dist,
    "entropy": compression.entropy_ncd,
    "fragile_ends_simil": edit.fragile_ends_simil,
    "jaccard": token.jaccard_dist,
    "jaro": edit.jaro_dist,
    "jaro_winkler": edit.jaro_winkler_dist,
    "levenshtein": edit.levenshtein_dist,
    "mmcwpa": edit.mmcwpa_dist,
    "ratcliff_obershelp": sequence.ratcliff_obershelp,
    "sorensen": token.sorensen_dist,
    "stemmatological_simil": edit.stemmatological_simil,
    "subseq_jaccard": token.subseq_jaccard_dist,
}

# TODO: accept other parameters via *kwargs
def distance(
    seqs: Sequence[Sequence[Hashable]],
    method: str = "levenshtein",
    normal: bool = False,
) -> float:
    """
    Computes the distance between sequences according to a specified method.

    This function acts as a wrapper to all the methods offered by the package,
    including those that are not properly "distances" but measures of
    "similarity" (that is, those that do not offer all the distance
    properties). It is intended as a single point of call for all the
    methods that are offered.

    Contrary to the individual methods that accept two sequence as arguments,
    this wrapper accepts a sequence of sequence, allowing to compute
    multiple distances.

    Examples
    *********

    .. code-block:: python

        >>> seqsim.distance(["abc", "bcde"])
        3.0
        >>> seqsim.distance(["abc", "bcde", "fgh"])
        3.3333333333333335

    :param seqs: A group of group of hashable elements to be compared.
        Currently, if more than two sequences are passed, it just returns
        the mean value of all pairwise comparisons, but this operation
        might change in the future at least for some methods.
    :param method: The method for comparison to be used. The list of
        methods, and the function they call, can be obtained from the
        keys of the `METHODS` dictionary exported by this module.
        Defaults to "levenshtein".
    :param normal: Whether to return a normalized score for the comparison
        in range [0..1]. Note that the function will accept a `True` value
        for all methods, but not all methods offer normalization and some
        method always return normalized values. In those cases, the
        standard value will be returned a warning message will be sent
        to the standard logger (which can be silenced as usual with the
        Python `logging` standard library. Defaults to `False`.
    :return: The distance score.
    """

    # Make sure we have at least two sequences
    if len(seqs) < 2:
        raise ValueError("At least two sequences are need for computation.")

    # Make sure the requested method is available
    if method not in METHODS:
        raise ValueError(f"Unknown or unsupported method `{method}.")

    # While we could use a combination also for two sequences, it takes a little
    # less code to check and do it directly, and it will make easier to
    # implement methods different than the mean in the future.
    if len(seqs) == 2:
        dist = METHODS[method](seqs[0], seqs[1], normal=normal)
    else:
        dists = [
            METHODS[method](seq_x, seq_y, normal=normal)
            for seq_x, seq_y in itertools.combinations(seqs, 2)
        ]
        dist = sum(dists) / len(seqs)

    return dist


# Build namespace
__all__ = [
    "distance",
    "METHODS",
]
