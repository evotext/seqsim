"""
seqsim __init__.py
"""

# Version of the seqsim package
__author__ = "Tiago Tresoldi"
__email__ = "tiago.tresoldi@lingfil.uu.se"
__version__ = "0.3"

# Import local modules
from . import distance
from . import similarity
from .common import set_seeds, collect_subseqs
from .ngrams import ngrams_iter, get_all_ngrams_by_order

# Build namespace
__all__ = ["set_seeds", "collect_subseqs", "ngrams_iter", "get_all_ngrams_by_order"]
