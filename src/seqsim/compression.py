"""
Module implementing various methods for similarity and distance from compression methods.

Most of the methods are commonly used in string comparison from normalized
compression distance, such as from Arithmetic Coding, but in this module we need
to make sure we can operate on arbitrary iterable data structures.
"""

# Import Python standard libraries
from typing import Hashable, Sequence
import logging

# Import 3rd-party libraries
import textdistance

# Import other modules
from . import common

# TODO: multiple sequences?
# TODO: consider splitting in similarity/distance
# TODO: have a normalization method based on seq_x*2 and seq_y*2
def arith_ncd(
    seq_x: Sequence[Hashable], seq_y: Sequence[Hashable], normal: bool = False
) -> float:
    """
    Computes a distance between two sequences based on Arithmetic Coding.

    See: https://en.wikipedia.org/wiki/Arithmetic_coding

    Example
    ********

    .. code-block:: python

        >>> seqsim.compression.arith_ncd("abc", "bcde")
        1.2222222222222223

    :param seq_x: The first sequence to be compared.
    :param seq_y: The second sequence to be compared.
    :param normal: Dummy parameter, see comment above.
    :return: The Arithmetic Coding NCD between the two sequences.
    """

    # As the method uses .find, we need strings
    seq_x, seq_y = common.equivalent_string(seq_x, seq_y)

    if normal:
        logging.warning("Arithmetic Coding NCD cannot be normalized in range [0..1].")

    return textdistance.arith_ncd(seq_x, seq_y)


# TODO: multiple sequences?
# TODO: have a normalization method based on seq_x*2 and seq_y*2
def entropy_ncd(
    seq_x: Sequence[Hashable], seq_y: Sequence[Hashable], normal: bool = False
) -> float:
    """
    Computes a distance between two sequences based on entropy.

    See: https://en.wikipedia.org/wiki/Entropy_(information_theory)

    Example
    ********

    .. code-block:: python

        >>> seqsim.compression.entropy_ncd("abc", "bcde")
        0.21698794996929216

    :param seq_x: The first sequence to be compared.
    :param seq_y: The second sequence to be compared.
    :param normal: Dummy parameter, see comment above.
    :return: The Entropy NCD between the two sequences.
    """

    # As the method uses .find, we need strings
    seq_x, seq_y = common.equivalent_string(seq_x, seq_y)

    if normal:
        logging.warning("Entropy NCD cannot be normalized in range [0..1].")

    return textdistance.entropy_ncd(seq_x, seq_y)
