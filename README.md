# seqsim

[![PyPI](https://img.shields.io/pypi/v/seqsim.svg)](https://pypi.org/project/seqsim)
[![CI](https://github.com/evotext/seqsim/actions/workflows/main.yml/badge.svg)](https://github.com/evotext/seqsim/actions/workflows/main.yml)
[![Documentation Status](https://readthedocs.org/projects/seqsim/badge/?version=latest)](https://seqsim.readthedocs.io/en/latest/?badge=latest)

Python library for computing measures of distance and similarity for sequences of hashable data types.

![scriptorium](https://raw.githubusercontent.com/evotext/seqsim/main/docs/scriptorium_small.jpg)

While developed as a general-purpose library, `seqsim` is mostly designed for usage
in research within the field of cultural evolution, and particularly of the
cultural evolution of textual traditions. Some methods act as a thin-wrapper
to the standard Python library; some implementations were ported from
[textdistance](https://github.com/life4/textdistance), and the library has no
third-party dependencies.

## Installation

In any standard Python environment, `seqsim` can be installed with:

```bash
$ pip install seqsim
```

## Usage

The library offers different methods to compare sequences of arbitrary hashable elements.
It is possible to mix sequence and element types.

Full documentation is offered at [ReadTheDocs](https://seqsim.readthedocs.io/en/latest/?badge=latest) and
code with almost complete coverage is offered in the
[tests](https://github.com/evotext/seqsim/tree/main/tests). For most common usages,
a wrapper `.distance()` function can be used.

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
>>> seqsim.distance(["abcdeXXXXXfghij", "abcdefghij"], "bulk_delete", max_del_len=5)
1.0
```

All functions take the two sequences as positional arguments; every other
parameter (such as `normal`, which requests a value in range [0..1]) must be
passed by name. With more than two sequences, `distance()` returns the mean of
all pairwise comparisons.

### Distances, dissimilarities, and similarities

Function names state the mathematical properties of each measure:

- **`_dist`**: a true distance (metric), with non-negativity, symmetry,
  identity of indiscernibles (`d(x, y) == 0` only if `x == y`), and the
  triangle inequality. For edit distances these hold for the raw values;
  normalized values do not satisfy the triangle inequality.
- **`_dissim`**: a dissimilarity, where identical sequences score `0.0` and
  higher values indicate more different sequences, but where the metric
  properties are not all guaranteed.
- **`_simil`**: a similarity, where higher values indicate more similar
  sequences.

All measures are symmetric. Comparing two empty sequences gives `0.0`, and
comparing an empty sequence with a non-empty one gives the maximum value
(`1.0` for measures in range [0..1]).

| Method (`distance()` key) | Function | Identity of indiscernibles | Triangle inequality | Range |
|---|---|---|---|---|
| `levenshtein` | `edit.levenshtein_dist` | yes | yes | [0..max length] |
| `damerau` | `edit.damerau_dist` | yes | yes | [0..max length] |
| `bulk_delete` | `edit.bulk_delete_dist` | yes | yes | [0..max length] |
| `osa` | `edit.osa_dissim` | yes | no | [0..max length] |
| `fragile_ends` | `edit.fragile_ends_dissim` | yes | no | [0..max length] |
| `stemmatological` | `edit.stemmatological_dissim` | yes | no | [0..max length] |
| `jaro` | `edit.jaro_dissim` | yes | no | [0..1] |
| `jaro_winkler` | `edit.jaro_winkler_dissim` | yes | no | [0..1] |
| `mmcwpa` | `edit.mmcwpa_dissim` | yes | no | [0..1] |
| `birnbaum` | `edit.birnbaum_dissim` | yes | no | [0..1] |
| `ratcliff_obershelp` | `sequence.ratcliff_obershelp_dissim` | yes | no | [0..1] |
| `subseq_jaccard` | `token.subseq_jaccard_dissim` | yes | not proven | [0..1] |
| `jaccard` | `token.jaccard_dissim` | no (ignores order and repetition) | yes | [0..1] |
| `sorensen` | `token.sorensen_dissim` | no (ignores order) | no | [0..1] |
| `entropy_ncd` | `compression.entropy_ncd_dissim` | no (ignores order) | no | [0..1] |
| `lzma_ncd` | `compression.lzma_ncd_dissim` | no (identical short sequences score above 0) | not guaranteed | [0..1] normalized |

`edit.birnbaum_simil()` is also available as a similarity score. The
"lzma_ncd" method, like any Normalized Compression Distance, is only meaningful
for sequences long enough to be compressed (dozens of elements or more).

## Demonstration

The table below, generated with `extra/readme_compare.py`, shows the results of
all methods for two pairs of sequences.

| Method             | Function                             |   "kitten" / "sitting" |   normalized |   (1, 2, 3, 4) / (3, 4, 2, 1) |   normalized |
|--------------------|--------------------------------------|------------------------|--------------|-------------------------------|--------------|
| birnbaum           | `edit.birnbaum_dissim`               |                 0.7500 |       0.7500 |                        0.7000 |       0.7000 |
| bulk_delete        | `edit.bulk_delete_dist`              |                 3.0000 |       0.4286 |                        2.0000 |       0.5000 |
| damerau            | `edit.damerau_dist`                  |                 3.0000 |       0.4286 |                        4.0000 |       1.0000 |
| entropy_ncd        | `compression.entropy_ncd_dissim`     |                 0.1013 |       0.1013 |                        0.0000 |       0.0000 |
| fragile_ends       | `edit.fragile_ends_dissim`           |                 2.5000 |       0.3571 |                        4.0000 |       1.0000 |
| jaccard            | `token.jaccard_dissim`               |                 0.5714 |       0.5714 |                        0.0000 |       0.0000 |
| jaro               | `edit.jaro_dissim`                   |                 0.2540 |       0.2540 |                        0.5000 |       0.5000 |
| jaro_winkler       | `edit.jaro_winkler_dissim`           |                 0.2540 |       0.2540 |                        0.5000 |       0.5000 |
| levenshtein        | `edit.levenshtein_dist`              |                 3.0000 |       0.4286 |                        4.0000 |       1.0000 |
| lzma_ncd           | `compression.lzma_ncd_dissim`        |                 0.6364 |       0.6364 |                        0.5000 |       0.5000 |
| mmcwpa             | `edit.mmcwpa_dissim`                 |                 0.5135 |       0.5135 |                        0.3876 |       0.3876 |
| osa                | `edit.osa_dissim`                    |                 3.0000 |       0.4286 |                        4.0000 |       1.0000 |
| ratcliff_obershelp | `sequence.ratcliff_obershelp_dissim` |                 0.3846 |       0.3846 |                        0.5000 |       0.5000 |
| sorensen           | `token.sorensen_dissim`              |                 0.3846 |       0.3846 |                        0.0000 |       0.0000 |
| stemmatological    | `edit.stemmatological_dissim`        |                 2.5000 |       0.3571 |                        2.0000 |       0.5000 |
| subseq_jaccard     | `token.subseq_jaccard_dissim`        |                 0.9549 |       0.9549 |                        0.8600 |       0.8600 |

## Changelog

See [CHANGELOG.md](CHANGELOG.md). Version 0.4.0 renames most functions and
changes several results; see the changelog for a migration guide.

## Community guidelines

While the authors can be contacted directly for support, it is recommended that third 
parties use GitHub standard features, such as issues and pull requests, to contribute, 
report problems, or seek support.

Contributing guidelines, including a code of conduct, can be found in the
`CONTRIBUTING.md` file.

## Authors and citation

The library is developed in the context of "[Cultural Evolution of Text](https://www.evotext.se)",
project, with funding from the Riksbankens Jubileumsfond (grant agreement ID:
[MXM19-1087:1](https://www.rj.se/en/anslag/2019/cultural-evolution-of-texts/)).

If you use `seqsim`, please cite it as:

> Tresoldi, Tiago; Maurits, Luke; Dunn, Michael. (2021). seqsim, a library
> for computing measures of distance and similarity for sequences of hashable data
> types. Version 0.4.0. Uppsala: Uppsala universitet.
> Available at: https://github.com/evotext/seqsim

In BibTeX:

```
@misc{Tresoldi2021seqsim,
  author = {Tresoldi, Tiago; Maurits, Luke; Dunn, Michael},
  title = {seqsim, a library for computing measures of distance and similarity for sequences of hashable data types. Version 0.4.0},
  howpublished = {\url{https://github.com/evotext/seqsim}},
  address = {Uppsala},
  publisher = {Uppsala universitet},
  year = {2021},
}
```

## References

The image at the top of this file is derived from Yves de Saint-Denis, *Vie et martyre de saint
Denis et de ses compagnons, versions latine et française*. It is available in high
resolution from [Bibliothèque nationale de France, Département des Manuscrits, Français 2090,
fol. 12v.](http://gallica.bnf.fr/ark:/12148/btv1b8447296x/f30.item)

References to the various implementation are available in the source code comments and in
the [online documentation](https://seqsim.readthedocs.io/en/latest/?badge=latest).
