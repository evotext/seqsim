"""
Module for collecting ngrams on sequences of arbitrary elements.

Most of the code in this module follows the original implementation by Tiago Tresoldi
for the `lingpy` library, later moved into the independent `lpngram` package.
"""

from itertools import chain
from typing import Sequence


# This method with zip, besides returning an iterator as desired, is faster
# than both the previous lingpy implementation and the one in NLTK; as this is
# the core of the ngram methods, it is important to have at least this
# primitive as fast as possible. This is intentionally not defaulting to any
# value for the order, so that users won't confuse a given order to all
# orders up to and including the given one.
# TODO: typing for pad_symbol and return
def ngrams_iter(sequence: Sequence, order: int, pad_symbol="$$$"):
    """
    Build an iterator for collecting all ngrams of a given order.

    The sequence can optionally be padded with boundary symbols which are
    equal for before and and after sequence boundaries.

    :param sequence: The sequence from which the ngrams will be collected.
    :param order: The order of the ngrams to be collected.
    :param pad_symbol: An optional symbol to be used as start-of- and end-of-sequence
        boundaries. The same symbol is used for both boundaries. Must be a
        value different from None, defaults to "$$$".
    """

    # Makes sure the sequence is a tuple, for faster computation, and pad it if needed.
    seq = tuple(sequence)
    if pad_symbol is not None:
        seq = chain((pad_symbol,) * (order - 1), seq, (pad_symbol,) * (order - 1))
        seq = tuple(seq)

    # We generate the collection of ngrams for counting occurrences with Python
    # `zip()` function, so we can rely on internal C-code for speeding things
    # up. What we do is to build a list of arguments for the function by using a
    # list comprehension over variable `order` and, then, decompose such list
    # when passing it to the function. The list comprehension is built so that,
    # for a sequence such as the characters in "Markov" (here without
    # boundaries, for simplicity) and an order of 3, we'll have the following elements,
    # from which we zip all possible combinations:
    #   [['M', 'a', 'r', 'k', 'o', 'v'],
    #    ['a', 'r', 'k', 'o', 'v'],
    #    ['r', 'k', 'o', 'v'],
    #    ['k', 'o', 'v']]
    yield from zip(*[seq[i:] for i in range(order)])
