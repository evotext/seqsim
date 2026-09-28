.. seqsim documentation master file, created by
   sphinx-quickstart on Wed Mar  3 11:12:32 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to seqsim's documentation!
==================================

|PyPI| |CI| |Documentation Status|

Python library for computing measures of distance and similarity for
sequences of hashable data types.

.. figure:: https://raw.githubusercontent.com/evotext/seqsim/main/docs/scriptorium_small.jpg
   :alt: scriptorium

While developed as a general-purpose library, ``seqsim`` is mostly
designed for usage in research within the field of cultural evolution,
and particularly of the cultural evolution of textual traditions. Some
methods act as a thin-wrapper to the standard Python library; some
implementations were ported from `textdistance`_, and the library has no
third-party dependencies.

Installation
------------

In any standard Python environment, ``seqsim`` can be installed with:

.. code:: bash

   $ pip install seqsim

Usage
-----

The library offers different methods to compare sequences of arbitrary
hashable elements. It is possible to mix sequence and element types.
For most common usages, a wrapper ``distance()`` function can be used.

.. code:: python

   >>> import seqsim
   >>> seqsim.edit.levenshtein_dist("kitten", "string")
   5.0
   >>> seqsim.edit.levenshtein_dist("kitten", "string", normal=True)
   0.8333333333333334
   >>> seqsim.edit.damerau_dist(["in", "the", "beginning"], ["the", "in", "beginning"])
   1.0
   >>> seqsim.distance(["kitten", "sitting", "fitting"], "jaro_winkler")
   0.20105820105820105
   >>> seqsim.distance(["abcdeXXXXXfghij", "abcdefghij"], "bulk_delete", max_del_len=5)
   1.0

All functions take the two sequences as positional arguments; every other
parameter (such as ``normal``, which requests a value in range [0..1]) must
be passed by name.

Function names state the mathematical properties of each measure:

- ``_dist``: a true distance (metric), with non-negativity, symmetry,
  identity of indiscernibles, and the triangle inequality (for edit
  distances, on the raw values);
- ``_dissim``: a dissimilarity, where identical sequences score ``0.0`` and
  higher values indicate more different sequences, but where the metric
  properties are not all guaranteed;
- ``_simil``: a similarity, where higher values indicate more similar
  sequences.

All measures are symmetric. The properties of each method, a comparison
table, and the changelog (including a migration guide for version 0.4.0)
are available in the `README`_ and `CHANGELOG`_.

Authors and citation
--------------------

The library is developed in the context of “`Cultural Evolution of
Text`_”, project, with funding from the Riksbankens Jubileumsfond (grant
agreement ID: `MXM19-1087:1`_).

If you use ``seqsim``, please cite it as:

   Tresoldi, Tiago; Maurits, Luke; Dunn, Michael. (2021). seqsim, a
   library for computing measures of distance and similarity for
   sequences of hashable data types. Version 0.4.0. Uppsala: Uppsala universitet.
   Available at: https://github.com/evotext/seqsim

In BibTeX:

::

   @misc{Tresoldi2021seqsim,
     author = {Tresoldi, Tiago; Maurits, Luke; Dunn, Michael},
     title = {seqsim, a library for computing measures of distance and similarity for sequences of hashable data types. Version 0.4.0},
     howpublished = {\url{https://github.com/evotext/seqsim}},
     address = {Uppsala},
     publisher = {Uppsala universitet},
     year = {2021},
   }

References
----------

The image at the top of this file is derived from Yves de Saint-Denis,
*Vie et martyre de saint Denis et de ses compagnons, versions latine et
française*. It is available in high resolution from `Bibliothèque
nationale de France, Département des Manuscrits, Français 2090, fol.
12v.`_

References to the various implementation are available in the source
code comments and in the `online documentation`_.

.. _Cultural Evolution of Text: https://www.evotext.se
.. _`MXM19-1087:1`: https://www.rj.se/en/anslag/2019/cultural-evolution-of-texts/
.. _Bibliothèque nationale de France, Département des Manuscrits, Français 2090, fol. 12v.: http://gallica.bnf.fr/ark:/12148/btv1b8447296x/f30.item
.. _online documentation: https://seqsim.readthedocs.io/en/latest/?badge=latest

.. _ReadTheDocs: https://seqsim.readthedocs.io/en/latest/?badge=latest
.. _tests: https://github.com/evotext/seqsim/tree/main/tests
.. _README: https://github.com/evotext/seqsim/blob/main/README.md
.. _CHANGELOG: https://github.com/evotext/seqsim/blob/main/CHANGELOG.md

.. _textdistance: https://github.com/life4/textdistance

.. |PyPI| image:: https://img.shields.io/pypi/v/seqsim.svg
   :target: https://pypi.org/project/seqsim
.. |CI| image:: https://github.com/evotext/seqsim/actions/workflows/main.yml/badge.svg
   :target: https://github.com/evotext/seqsim/actions/workflows/main.yml
.. |Documentation Status| image:: https://readthedocs.org/projects/seqsim/badge/?version=latest
   :target: https://seqsim.readthedocs.io/en/latest/?badge=latest

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   Modules <source/modules.rst>


Indices and tables
==================

* :ref:`genindex`
* :ref:`search`
