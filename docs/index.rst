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
methods act as a thin-wrapper to either the standard Python library or
of to other libraries such as `textdistance`_.

Installation
------------

In any standard Python environment, ``seqsim`` can be installed with:

.. code:: bash

   $ pip install seqsim

Usage
-----

The library offers different methods to compare sequences of arbitrary
hashable elements. It is possible to mix sequence and element types.

Full documentation is offered at `ReadTheDocs`_ and code with almost
complete coverage is offered in the `tests`_. For most common usages, a
wrapper ``.distance()`` function can be used.

.. code:: python

   >>> import seqsim
   >>> seqsim.edit.levenshtein_dist("kitten", "string")
   5
   >>> seqsim.edit.levenshtein_dist("kitten", "string", normal=True)
   >>> 0.8333333333333334
   >>> seqsim.sequence.ratcliff_obershelp([1,2,3,4], [2,4,3,5])
   0.5
   >>> seqsim.compression.entropy_ncd([1,2,3,4], [2,4,3,5])
   0.08333333333333333

.. _ReadTheDocs: https://seqsim.readthedocs.io/en/latest/?badge=latest
.. _tests: https://github.com/evotext/seqsim/tree/main/tests

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
