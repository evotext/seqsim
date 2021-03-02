"""
seqsim __init__.py

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

# Import local modules
from . import edit
from . import token
from .common import set_seeds, collect_subseqs
from .ngrams import ngrams_iter, get_all_ngrams_by_order

# Build namespace
__all__ = ["set_seeds", "collect_subseqs", "ngrams_iter", "get_all_ngrams_by_order"]
