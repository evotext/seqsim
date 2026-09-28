"""
Module implementing methods for sequence dissimilarity based on compression.

The methods in this module compute a Normalized Compression Distance (NCD),
defined by Cilibrasi & Vitányi (2005) as

    NCD(x, y) = (C(xy) - min(C(x), C(y))) / max(C(x), C(y))

where `C(s)` is the size of `s` after compression. Each method differs in the
"compressor" it uses. All methods operate on sequences of arbitrary hashable
elements and are implemented in this package, so their results do not depend
on third-party libraries.

References
**********

Cilibrasi, Rudi; Vitányi, Paul M.B. (2005). "Clustering by compression". IEEE
Transactions on Information Theory 51 (4): 1523–1545.
"""

# Import Python standard libraries
from collections import Counter
from typing import Callable, Hashable, Sequence, Tuple
import lzma
import math

# Minimum and maximum dictionary sizes for LZMA compression; the dictionary
# only needs to cover the data being compressed, and allocating the default
# (64 MiB for preset 9) would make each call orders of magnitude slower
_LZMA_MIN_DICT = 1 << 16
_LZMA_MAX_DICT = 1 << 30


def _encode_bytes(
    seq_x: Sequence[Hashable], seq_y: Sequence[Hashable]
) -> Tuple[bytes, bytes]:
    """
    Maps a pair of sequences to equivalent byte strings.

    Each distinct element is mapped to a fixed-width code (one byte if there
    are at most 256 distinct elements, two bytes if at most 65536, and so on),
    so that general-purpose compressors can operate on sequences of any
    hashable elements. The mapping is shared by both sequences and does not
    depend on the order of the arguments.
    """

    elements = sorted(
        dict.fromkeys([*seq_x, *seq_y]),
        key=lambda element: (type(element).__qualname__, repr(element)),
    )
    width = max(1, math.ceil(math.log(max(len(elements), 1), 256)))
    codes = {
        element: idx.to_bytes(width, "big") for idx, element in enumerate(elements)
    }

    return (
        b"".join(codes[element] for element in seq_x),
        b"".join(codes[element] for element in seq_y),
    )


def _ncd(
    comp_x: float, comp_y: float, comp_xy: float, comp_yx: float
) -> float:
    """
    Computes the NCD from the compressed sizes.

    The size of the concatenation is the smallest of both orders, so that the
    result is symmetric.
    """

    max_comp = max(comp_x, comp_y)
    if max_comp == 0:
        return 0.0

    return (min(comp_xy, comp_yx) - min(comp_x, comp_y)) / max_comp


def _lzma_size(data: bytes) -> int:
    """
    Returns the size of the raw LZMA2 compression of `data`.
    """

    # Using the raw format (without the .xz container) avoids headers and
    # checksums that would dominate the compressed size of short sequences
    dict_size = min(max(len(data), _LZMA_MIN_DICT), _LZMA_MAX_DICT)
    filters = [
        {
            "id": lzma.FILTER_LZMA2,
            "preset": 9 | lzma.PRESET_EXTREME,
            "dict_size": dict_size,
        }
    ]

    return len(lzma.compress(data, format=lzma.FORMAT_RAW, filters=filters))


def lzma_ncd(
    seq_x: Sequence[Hashable], seq_y: Sequence[Hashable], normal: bool = False
) -> float:
    """
    Computes the Normalized Compression Distance using the LZMA compressor.

    Sequences are mapped to byte strings (one fixed-width code per distinct
    element) and compressed with a raw LZMA2 stream from the Python standard
    library. As with any NCD, the result is only meaningful for sequences
    long enough for the compressor to exploit repetition (in the order of
    dozens of elements or more): identical short sequences do not score zero.

    The raw NCD can be slightly larger than 1.0 due to compressor overhead;
    when `normal` is set, the value is clipped to the range [0..1].

    Example
    ********

    .. code-block:: python

        >>> seqsim.compression.lzma_ncd("abc", "bcde")
        0.5

    :param seq_x: The first sequence to be compared.
    :param seq_y: The second sequence to be compared.
    :param normal: Whether to clip the result to the range [0..1].
    :return: The LZMA NCD between the two sequences.
    """

    bytes_x, bytes_y = _encode_bytes(seq_x, seq_y)
    ncd = _ncd(
        _lzma_size(bytes_x),
        _lzma_size(bytes_y),
        _lzma_size(bytes_x + bytes_y),
        _lzma_size(bytes_y + bytes_x),
    )

    if normal:
        return min(max(ncd, 0.0), 1.0)

    return ncd


def _entropy_size(seq: Sequence[Hashable]) -> float:
    """
    Returns the "compressed size" used by the entropy NCD.

    The size is one plus the Shannon entropy (in bits per symbol) of the
    distribution of elements in the sequence, following the definition in the
    `textdistance` library (version 4.5), from which this method was ported.
    """

    total = len(seq)
    entropy = 0.0
    if total:
        for count in Counter(seq).values():
            prob = count / total
            entropy -= prob * math.log2(prob)

    return 1.0 + entropy


def entropy_ncd(
    seq_x: Sequence[Hashable], seq_y: Sequence[Hashable], normal: bool = False
) -> float:
    """
    Computes a Normalized Compression Distance based on entropy.

    The "compressed size" of a sequence is one plus the Shannon entropy of the
    distribution of its elements. As such, the method only considers the
    frequency of elements, not their order: any two sequences with the same
    element frequencies (e.g., `"ab"` and `"ba"`, or `"a"` and `"aaaa"`) have a
    distance of zero. The results are always in range [0..1], so `normal` has
    no effect.

    This is a port of the `EntropyNCD` method of the `textdistance` library.

    Example
    ********

    .. code-block:: python

        >>> seqsim.compression.entropy_ncd("abc", "bcde")
        0.21698794996929216

    References
    ***********

    Shannon, C.E., Weaver, W. (1949) The Mathematical Theory of Communication, Univ of
    Illinois Press. ISBN 0-252-72548-4

    :param seq_x: The first sequence to be compared.
    :param seq_y: The second sequence to be compared.
    :param normal: Ignored, as results are always in range [0..1].
    :return: The Entropy NCD between the two sequences.
    """

    concat = [*seq_x, *seq_y]
    size_xy = _entropy_size(concat)

    return _ncd(_entropy_size(seq_x), _entropy_size(seq_y), size_xy, size_xy)
