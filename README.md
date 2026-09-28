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

The [documentation](https://seqsim.readthedocs.io) includes a guide to the
concepts and to choosing a method, tutorials on comparing witnesses of a text
and the contents of manuscripts, a case study on a real textual tradition,
and a detailed description of every measure. For most common usages, a
wrapper `.distance()` function can be used.

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

All methods of `distance()` are symmetric. Comparing two empty sequences gives `0.0`, and
comparing an empty sequence with a non-empty one gives the maximum value
(`1.0` for measures in range [0..1]).

| Method (`distance()` key) | Function | Identity of indiscernibles | Triangle inequality | Range |
|---|---|---|---|---|
| `levenshtein` | `edit.levenshtein_dist` | yes | yes | [0..max length] |
| `damerau` | `edit.damerau_dist` | yes | yes | [0..max length] |
| `indel` | `edit.indel_dist` | yes | yes | [0..len(x) + len(y)] |
| `bulk_delete` | `edit.bulk_delete_dist` | yes | yes | [0..max length] |
| `lcs` | `edit.lcs_dist` | yes | yes | [0..1] |
| `levenshtein_gld` | `edit.levenshtein_gld_dist` | yes | yes | [0..1] |
| `damerau_gld` | `edit.damerau_gld_dist` | yes | yes | [0..1] |
| `indel_gld` | `edit.indel_gld_dist` | yes | yes | [0..1] |
| `levenshtein_ned` | `edit.levenshtein_ned_dist` | yes | yes | [0..1] |
| `ulam` | `order.ulam_dist` | yes | yes | [0..len(x) + len(y)] |
| `osa` | `edit.osa_dissim` | yes | no | [0..max length] |
| `fragile_ends` | `edit.fragile_ends_dissim` | yes | no | [0..max length] |
| `stemmatological` | `edit.stemmatological_dissim` | yes | no | [0..max length] |
| `nw` | `alignment.nw_dissim` | yes | with metric costs and no gap opening | [0..gap costs] |
| `block_move` | `edit.block_move_dissim` | yes | no | [0..max length + 1] |
| `kendall_tau` | `order.kendall_tau_dissim` | yes | no (near metric) | [0..number of pairs] |
| `footrule` | `order.footrule_dissim` | yes | with a fixed `ell` | [0..] |
| `cayley` | `order.cayley_dissim` | yes | no | [0..len(x) + len(y)] |
| `block_interchange` | `order.block_interchange_dissim` | yes | not proven | [0..len(x) + len(y)] |
| `jaro` | `edit.jaro_dissim` | yes | no | [0..1] |
| `jaro_winkler` | `edit.jaro_winkler_dissim` | yes | no | [0..1] |
| `mmcwpa` | `edit.mmcwpa_dissim` | yes | no | [0..1] |
| `birnbaum` | `edit.birnbaum_dissim` | yes | no | [0..1] |
| `ratcliff_obershelp` | `sequence.ratcliff_obershelp_dissim` | yes | no | [0..1] |
| `subseq_jaccard` | `token.subseq_jaccard_dissim` | yes | not proven | [0..1] |
| `gst` | `edit.gst_dissim` | no (ignores the order of tiles) | no | [0..1] |
| `breakpoint` | `order.breakpoint_dissim` | no (with repeated elements) | yes | [0..] |
| `qgram` | `token.qgram_dissim` | no | yes | [0..] |
| `jaccard` | `token.jaccard_dissim` | no (ignores order and repetition) | yes | [0..1] |
| `sorensen` | `token.sorensen_dissim` | no (ignores order) | no | [0..1] |
| `entropy_ncd` | `compression.entropy_ncd_dissim` | no (ignores order) | no | [0..1] |
| `lzma_ncd` | `compression.lzma_ncd_dissim` | no (identical short sequences score above 0) | not guaranteed | [0..1] normalized |
| `lz76` | `compression.lz76_dissim` | no (e.g. `aa` and `aaa`) | no | [0..1] normalized |

All methods accept `normal=True`, which returns a value in range [0..1].

Other functions, outside `distance()`:

- similarities: `edit.birnbaum_simil`, `alignment.sw_simil` (local alignment),
  `alignment.monge_elkan_simil` (sequences of sequences, such as lists of
  titles), and `token.tversky_simil` (directional unless `alpha == beta`);
- `token.containment`, how much of one sequence is found in another
  (directional);
- `order.iebp_estimate`, an estimate of the number of transpositions
  separating the shared items of two sequences (Spencer et al., 2003).

Measures in the `order` module compare the order of shared items, such as the
texts in manuscripts with overlapping contents; repeated elements are matched
by occurrence. The "lzma_ncd" and "lz76" methods, like any compression-based measure, are only meaningful
for sequences long enough to be compressed (dozens of elements or more).

## Demonstration

The table below, generated with `extra/readme_compare.py`, shows the results of
all methods for two pairs of sequences.

| Method             | Function                             |   "kitten" / "sitting" |   normalized |   (1, 2, 3, 4) / (3, 4, 2, 1) |   normalized |
|--------------------|--------------------------------------|------------------------|--------------|-------------------------------|--------------|
| birnbaum           | `edit.birnbaum_dissim`               |                 0.7500 |       0.7500 |                        0.7000 |       0.7000 |
| block_interchange  | `order.block_interchange_dissim`     |                 5.0000 |       0.5556 |                        1.0000 |       0.2500 |
| block_move         | `edit.block_move_dissim`             |                 6.0000 |       0.7500 |                        4.0000 |       0.8000 |
| breakpoint         | `order.breakpoint_dissim`            |                 5.5000 |       0.7333 |                        4.0000 |       0.8000 |
| bulk_delete        | `edit.bulk_delete_dist`              |                 3.0000 |       0.4286 |                        2.0000 |       0.5000 |
| cayley             | `order.cayley_dissim`                |                 5.0000 |       0.5556 |                        3.0000 |       0.7500 |
| damerau            | `edit.damerau_dist`                  |                 3.0000 |       0.4286 |                        4.0000 |       1.0000 |
| damerau_gld        | `edit.damerau_gld_dist`              |                 0.3750 |       0.3750 |                        0.6667 |       0.6667 |
| entropy_ncd        | `compression.entropy_ncd_dissim`     |                 0.1013 |       0.1013 |                        0.0000 |       0.0000 |
| footrule           | `order.footrule_dissim`              |                21.0000 |       0.3818 |                        8.0000 |       0.4000 |
| fragile_ends       | `edit.fragile_ends_dissim`           |                 2.5000 |       0.3571 |                        4.0000 |       1.0000 |
| gst                | `edit.gst_dissim`                    |                 0.5385 |       0.5385 |                        0.5000 |       0.5000 |
| indel              | `edit.indel_dist`                    |                 5.0000 |       0.3846 |                        4.0000 |       0.5000 |
| indel_gld          | `edit.indel_gld_dist`                |                 0.5556 |       0.5556 |                        0.6667 |       0.6667 |
| jaccard            | `token.jaccard_dissim`               |                 0.5714 |       0.5714 |                        0.0000 |       0.0000 |
| jaro               | `edit.jaro_dissim`                   |                 0.2540 |       0.2540 |                        0.5000 |       0.5000 |
| jaro_winkler       | `edit.jaro_winkler_dissim`           |                 0.2540 |       0.2540 |                        0.5000 |       0.5000 |
| kendall_tau        | `order.kendall_tau_dissim`           |                23.0000 |       0.5111 |                        5.0000 |       0.5000 |
| lcs                | `edit.lcs_dist`                      |                 0.4286 |       0.4286 |                        0.5000 |       0.5000 |
| levenshtein        | `edit.levenshtein_dist`              |                 3.0000 |       0.4286 |                        4.0000 |       1.0000 |
| levenshtein_gld    | `edit.levenshtein_gld_dist`          |                 0.3750 |       0.3750 |                        0.6667 |       0.6667 |
| levenshtein_ned    | `edit.levenshtein_ned_dist`          |                 0.4286 |       0.4286 |                        0.6667 |       0.6667 |
| lz76               | `compression.lz76_dissim`            |                 0.5000 |       0.5000 |                        0.5000 |       0.5000 |
| lzma_ncd           | `compression.lzma_ncd_dissim`        |                 0.6364 |       0.6364 |                        0.5000 |       0.5000 |
| mmcwpa             | `edit.mmcwpa_dissim`                 |                 0.5135 |       0.5135 |                        0.3876 |       0.3876 |
| nw                 | `alignment.nw_dissim`                |                 3.0000 |       0.2308 |                        4.0000 |       0.5000 |
| osa                | `edit.osa_dissim`                    |                 3.0000 |       0.4286 |                        4.0000 |       1.0000 |
| qgram              | `token.qgram_dissim`                 |                11.0000 |       0.7333 |                        8.0000 |       0.8000 |
| ratcliff_obershelp | `sequence.ratcliff_obershelp_dissim` |                 0.3846 |       0.3846 |                        0.5000 |       0.5000 |
| sorensen           | `token.sorensen_dissim`              |                 0.3846 |       0.3846 |                        0.0000 |       0.0000 |
| stemmatological    | `edit.stemmatological_dissim`        |                 2.5000 |       0.3571 |                        2.0000 |       0.5000 |
| subseq_jaccard     | `token.subseq_jaccard_dissim`        |                 0.9549 |       0.9549 |                        0.8600 |       0.8600 |
| ulam               | `order.ulam_dist`                    |                 5.0000 |       0.5556 |                        2.0000 |       0.5000 |

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
