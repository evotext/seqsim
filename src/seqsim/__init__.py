"""
Main module of the `seqsim` package.

We follow the mathematical definitions for distinguishing between measures of
"distance", "dissimilarity", and "similarity", stated in the name of each
function:

  * `_dist`: a true distance (metric), with the properties of
      * non-negativity: d(x,y) >= 0
      * symmetry: d(x,y) = d(y,x)
      * identity of indiscernibles: d(x,y) = 0 <=> x = y
      * triangle inequality: d(x,z) <= d(x,y) + d(y,z)
  * `_dissim`: a dissimilarity, where 0.0 indicates identical sequences
    and higher values more different ones, but for which the properties
    above are not all guaranteed;
  * `_simil`: a similarity, where higher values indicate more similar
    sequences.
"""

# Version of the `seqsim` package
__author__ = "Tiago Tresoldi, Luke Maurits, Michael Dunn"
__email__ = "tiago.tresoldi@lingfil.uu.se"
__version__ = "0.4.0"

# Import Python standard libraries
from collections.abc import Sequence as SequenceABC
from typing import Any, Hashable, Iterable
import itertools

# Import local modules
from . import alignment
from . import edit
from . import order
from . import token
from . import sequence
from . import compression
from . import tradition
from .ngrams import ngrams_iter, get_all_ngrams_by_order
from ._measure import MeasureInfo, measures, methods

# All methods available through `distance()`, mapped to their functions and
# sorted by name. The mapping is built from the declarations of the measures
# (see `measures()`), and is a convenient single point of reference for
# users. All methods return 0.0 for identical sequences, with higher values
# for more different sequences.
METHODS = methods()


def distance(
    seqs: Iterable[Iterable[Hashable]],
    method: str = "levenshtein",
    *,
    normal: bool = False,
    **kwargs: Any,
) -> float:
    """
    Computes the distance between sequences according to a specified method.

    This function acts as a wrapper to all the methods offered by the package,
    including those that are not properly "distances" but dissimilarities
    (that is, those that do not offer all the distance properties). It is
    intended as a single point of call for all the methods that are offered.

    Contrary to the individual methods that accept two sequence as arguments,
    this wrapper accepts a collection of sequences, allowing to compute
    multiple distances.

    Examples
    *********

    .. code-block:: python

        >>> seqsim.distance(["abc", "bcde"])
        3.0
        >>> seqsim.distance(["abc", "bcde", "fgh"])
        3.3333333333333335
        >>> seqsim.distance(["abcdeXXXXXfghij", "abcdefghij"], "bulk_delete", max_del_len=5)
        1.0

    :param seqs: A collection of at least two sequences of hashable elements to
        be compared. Any iterable is accepted (including generators), but a
        single string is rejected, as it is most likely a mistake. If more than
        two sequences are passed, the mean of all pairwise comparisons is
        returned, but this operation might change in the future at least for
        some methods.
    :param method: The method for comparison to be used. The list of
        methods, and the function they call, can be obtained from the
        keys of the `METHODS` dictionary exported by this module.
        Defaults to "levenshtein".
    :param normal: Whether to return a normalized score for the comparison
        in range [0..1]. All methods accept this parameter; for methods whose
        results are always in range [0..1] it has no effect. Defaults to
        `False`.
    :param kwargs: Additional keyword arguments passed to the method, such as
        `max_del_len` for "bulk_delete".
    :return: The distance score.
    """

    # Make sure the requested method is available
    if method not in METHODS:
        raise ValueError(f"Unknown or unsupported method `{method}`.")

    # Reject single strings, which would be interpreted as a collection of
    # single-character sequences
    if isinstance(seqs, (str, bytes)):
        raise TypeError(
            "`seqs` must be a collection of sequences, not a single string; "
            "use e.g. `distance([seq_x, seq_y])`."
        )

    # Accept any iterable of iterables, making sure we have sequences
    seqs = [seq if isinstance(seq, SequenceABC) else tuple(seq) for seq in seqs]

    # Make sure we have at least two sequences
    if len(seqs) < 2:
        raise ValueError("At least two sequences are need for computation.")

    func = METHODS[method]
    dists = [
        func(seq_x, seq_y, normal=normal, **kwargs)
        for seq_x, seq_y in itertools.combinations(seqs, 2)
    ]

    return float(sum(dists) / len(dists))


# Build namespace
__all__ = [
    "distance",
    "METHODS",
    "measures",
    "MeasureInfo",
    "alignment",
    "edit",
    "order",
    "token",
    "sequence",
    "compression",
    "tradition",
    "ngrams_iter",
    "get_all_ngrams_by_order",
]
