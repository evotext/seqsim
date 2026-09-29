# seqsim

[![PyPI](https://img.shields.io/pypi/v/seqsim.svg)](https://pypi.org/project/seqsim)
[![Python](https://img.shields.io/pypi/pyversions/seqsim.svg)](https://pypi.org/project/seqsim)
[![License](https://img.shields.io/pypi/l/seqsim.svg)](https://github.com/tresoldi/seqsim/blob/main/LICENSE)
[![CI](https://github.com/tresoldi/seqsim/actions/workflows/main.yml/badge.svg)](https://github.com/tresoldi/seqsim/actions/workflows/main.yml)
[![Documentation Status](https://readthedocs.org/projects/seqsim/badge/?version=latest)](https://seqsim.readthedocs.io/en/latest/?badge=latest)
[![Code style: ruff](https://img.shields.io/badge/code%20style-ruff-261230.svg)](https://github.com/astral-sh/ruff)

![scriptorium](https://raw.githubusercontent.com/tresoldi/seqsim/main/docs/scriptorium_small.jpg)

`seqsim` is a Python library of measures of distance and similarity for
sequences of arbitrary hashable elements: characters, words, or, for
manuscripts, the texts they contain. It was designed for research on the
cultural evolution of textual traditions, but it is a general-purpose library.

- **A wide choice of measures.** Edit distances, alignment scores, token-based
  and order-based measures, and compression-based measures, in a single
  consistent interface.
- **Explicit mathematical properties.** Names state what a measure returns
  (`_dist`, `_dissim`, or `_simil`), and each claim (identity of
  indiscernibles, triangle inequality) is tested, with counterexamples for
  the claims that fail.
- **For sequences of anything.** Elements can be any hashable object, and
  different sequence types can be mixed.
- **Tools for whole traditions.** The `seqsim.tradition` subpackage prepares
  collections of manuscripts for phylogenetic software.
- **No dependencies.** Pure Python (3.10 or later), fully type-annotated.

## Installation

```bash
pip install seqsim
```

## Quick start

```python
>>> import seqsim
>>> seqsim.edit.levenshtein_dist("kitten", "string")
5.0
>>> seqsim.edit.levenshtein_dist("kitten", "string", normal=True)
0.8333333333333334
>>> seqsim.edit.damerau_dist(["in", "the", "beginning"], ["the", "in", "beginning"])
1.0
>>> seqsim.sequence.ratcliff_obershelp_dissim([1, 2, 3, 4], [2, 4, 3, 5])
0.5
>>> seqsim.distance(["kitten", "sitting", "fitting"], "jaro_winkler")
0.20105820105820105
```

All functions take the two sequences as positional arguments; every other
parameter, such as `normal=True` (which requests a value in range [0..1]),
must be passed by name. The `distance()` function is a single entry point to
the main measures by name, and with more than two sequences returns the mean
of all pairwise comparisons.

## Measures

| Module | Measures |
|---|---|
| `seqsim.edit` | Levenshtein, Damerau-Levenshtein, OSA, indel and LCS, normalized variants, Jaro and Jaro-Winkler, block moves, Greedy String Tiling, and measures for stemmatology |
| `seqsim.order` | Ulam, Kendall tau, Spearman footrule, Cayley, block interchange, breakpoints, and the IEBP estimate, for comparing the order of shared items |
| `seqsim.alignment` | Needleman-Wunsch, Smith-Waterman, and Monge-Elkan |
| `seqsim.token` | Jaccard, Sørensen-Dice, q-grams, Tversky, and containment |
| `seqsim.sequence` | Ratcliff-Obershelp |
| `seqsim.compression` | Normalized compression distances and Lempel-Ziv |

Function names state the mathematical properties of each measure:

- **`_dist`**: a true distance (metric), with non-negativity, symmetry,
  identity of indiscernibles, and the triangle inequality. For edit distances
  these hold for the raw values; normalized values do not satisfy the
  triangle inequality.
- **`_dissim`**: a dissimilarity, where identical sequences score `0.0` and
  higher values indicate more different sequences, but where the metric
  properties are not all guaranteed.
- **`_simil`**: a similarity, where higher values indicate more similar
  sequences.

The [documentation](https://seqsim.readthedocs.io) has the table of all
measures with their properties, a [guide to choosing a
method](https://seqsim.readthedocs.io/en/latest/choosing.html), and a
detailed description of each one. Measures are also listed programmatically by
`seqsim.measures()`.

## Traditions

The `seqsim.tradition` subpackage works on collections of witnesses that
share a reference frame, such as manuscripts of a collection of sayings. It
builds a consensus order of all items, decides where the absence of an item
from a witness is evidence and where it is missing data (lost leaves,
sections the witness never had, excerpts), codes content and adjacency
(order) as binary characters, and writes PHYLIP and NEXUS files for programs
such as IQ-TREE and MrBayes. It does not build trees.

```python
>>> from seqsim import tradition
>>> witnesses = {"W1": "abcdef", "W2": "abcdef", "W3": "abdcf", "W4": "abdcef"}
>>> frame = tradition.consensus_order(witnesses.values())
>>> label_of = dict.fromkeys(frame, "I")
>>> coverages = {w: tradition.coverage(seq, frame, label_of) for w, seq in witnesses.items()}
>>> content = tradition.content_characters(coverages, min_each=1)
>>> adjacency = tradition.adjacency_characters(witnesses, frame)
>>> print(tradition.to_phylip(tradition.concat_characters([content, adjacency])), end="")
4 5
W1 11010
W2 11010
W3 00101
W4 10101
```

## Documentation

Full documentation is at <https://seqsim.readthedocs.io>: a
[getting started](https://seqsim.readthedocs.io/en/latest/getting-started.html)
guide, tutorials on comparing witnesses and manuscript contents, a case study
on a real textual tradition, and the API reference.

## Contributing

Bug reports, questions, and pull requests are welcome on
[GitHub](https://github.com/tresoldi/seqsim/issues); see
[CONTRIBUTING.md](CONTRIBUTING.md) for the development setup and the
conventions for adding a method, and the [changelog](CHANGELOG.md) for the
history of releases. Version 0.4.0 renames most functions and changes several
results; the changelog has a migration guide.

## Citation

If you use `seqsim`, please cite it as:

> Tresoldi, Tiago (2026). *seqsim: measures of distance and similarity for
> sequences of hashable data types*. Version 0.4.0.
> <https://github.com/tresoldi/seqsim>

In BibTeX:

```bibtex
@software{tresoldi2026seqsim,
  author  = {Tresoldi, Tiago},
  title   = {seqsim: measures of distance and similarity for sequences of hashable data types},
  version = {0.4.0},
  year    = {2026},
  url     = {https://github.com/tresoldi/seqsim},
}
```

A machine-readable version is in [CITATION.cff](CITATION.cff), which GitHub
uses for its "Cite this repository" button. Please also cite the original
publications of the methods you use, which are listed in the documentation.

## Acknowledgements

The library was developed in the context of the project "[Cultural Evolution
of Texts](https://www.evotext.se)", with funding from the Riksbankens
Jubileumsfond (grant agreement ID
[MXM19-1087:1](https://www.rj.se/en/anslag/2019/cultural-evolution-of-texts/)).
Some implementations were ported from
[textdistance](https://github.com/life4/textdistance).
The "bulk delete", "fragile ends" and "stemmatological" measures were
developed by Luke Maurits and Michael Dunn for [Göransson et al.
(2023)](https://doi.org/10.1093/llc/fqac025).

The image at the top of this file is derived from Yves de Saint-Denis, *Vie et
martyre de saint Denis et de ses compagnons, versions latine et française*,
Bibliothèque nationale de France, Département des Manuscrits, Français 2090,
fol. 12v, available in high resolution from
[Gallica](http://gallica.bnf.fr/ark:/12148/btv1b8447296x/f30.item).

## License

MIT; see [LICENSE](LICENSE).
