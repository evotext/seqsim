# seqsim

**Measures of distance and similarity for sequences of any kind of element.**

`seqsim` compares sequences: the words of two witnesses of a text, the
characters of two spellings of a name, the list of texts copied in two
manuscripts, or any other ordered collection of Python objects. It offers
more than forty measures, from the classic Levenshtein distance to methods
developed for textual scholarship, such as the "stemmatological"
dissimilarity, and methods from the study of genome rearrangements that can
be applied to the order of texts in manuscripts.

While it is a general-purpose library, `seqsim` is designed with
stemmatology and the digital humanities in mind:

- **Any element, not only characters.** Sequences can hold words, lemmata,
  text identifiers, tuples, or any hashable Python object, and can mix types.
- **Honest names.** Every function states its mathematical properties:
  `_dist` for true distances (metrics), `_dissim` for dissimilarities, and
  `_simil` for similarities. This matters when the results feed
  distance-based methods for building trees and networks.
- **Measures for manuscript traditions.** Beyond edit distances, the
  library measures how the order of shared texts differs, how many blocks of
  texts were moved, how much of one collection is contained in another, and
  more.
- **No dependencies.** Pure Python, verified against reference
  implementations and published examples.

```python
>>> import seqsim
>>> witness_a = "in principio erat verbum".split()
>>> witness_b = "in principio erat sermo".split()
>>> seqsim.edit.levenshtein_dist(witness_a, witness_b)
1.0
>>> manuscript_1 = ["Vita Antonii", "Vita Pauli", "Vita Hilarionis", "Vita Malchi"]
>>> manuscript_2 = ["Vita Pauli", "Vita Hilarionis", "Vita Malchi", "Vita Antonii"]
>>> seqsim.order.ulam_dist(manuscript_1, manuscript_2)
1.0
```

## Where to start

- New to the library? Start with [Getting started](getting-started.md), then
  read [Concepts](concepts.md) for the ideas behind the measures.
- Looking for the right measure for a question? See
  [Choosing a method](choosing.md).
- Learn by example with the tutorials on
  [comparing witnesses of a text](tutorials/witnesses.md),
  [comparing the contents of manuscripts](tutorials/contents.md), and
  [identifying texts](tutorials/identifying.md), and with a
  [case study on a real tradition](tutorials/case-study.md).
- The [method pages](methods/index.md) describe every measure, with its
  properties and references.

```{toctree}
:hidden:
:maxdepth: 2
:caption: Guide

getting-started
concepts
choosing
```

```{toctree}
:hidden:
:maxdepth: 2
:caption: Tutorials

tutorials/witnesses
tutorials/contents
tutorials/identifying
tutorials/case-study
```

```{toctree}
:hidden:
:maxdepth: 2
:caption: Methods

methods/index
methods/edit
methods/order
methods/alignment
methods/token
methods/sequence
methods/compression
```

```{toctree}
:hidden:
:maxdepth: 1
:caption: Reference

api
references
changelog
contributing
citing
```
