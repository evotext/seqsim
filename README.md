# seqsim

[![PyPI](https://img.shields.io/pypi/v/seqsim.svg)](https://pypi.org/project/seqsim)
[![CI](https://github.com/evotext/seqsim/actions/workflows/main.yml/badge.svg)](https://github.com/evotext/seqsim/actions/workflows/main.yml)
[![Documentation Status](https://readthedocs.org/projects/seqsim/badge/?version=latest)](https://seqsim.readthedocs.io/en/latest/?badge=latest)

Python library for computing measures of distance and similarity for sequences of hashable data types.

![scriptorium](https://raw.githubusercontent.com/evotext/seqsim/main/docs/scriptorium_small.jpg)

While developed as a general-purpose library, `seqsim` is mostly designed for usage
in research within the field of cultural evolution, and particularly of the
cultural evolution of textual traditions. Some methods act as a thin-wrapper
to either the standard Python library or of to other libraries such as
[textdistance](https://github.com/life4/textdistance). 

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
5
>>> seqsim.edit.levenshtein_dist("kitten", "string", normal=True)
>>> 0.8333333333333334
>>> seqsim.sequence.ratcliff_obershelp([1,2,3,4], [2,4,3,5])
0.5
>>> seqsim.compression.entropy_ncd([1,2,3,4], [2,4,3,5])
0.08333333333333333
```

## Demonstration

The core of the library are the metrics for sequence distance/similarity on
arbitrary data types, as in the table below.

| Method                       |   "kitten" / "sitting" |   (1, 2, 3, 4) / (3, 4, 2, 1) |
|------------------------------|------------------------|-------------------------------|
| arith_ncd                    |               1.25     |                      0.888889 |
| arith_ncd_normal             |               1.25     |                      0.888889 |
| birnbaum_simil               |              10        |                      5        |
| birnbaum_simil_normal        |               0.3125   |                      0.5      |
| birnbaum                     |               0.565217 |                      0.5      |
| birnbaum_normal              |               0.565217 |                      0.5      |
| bulk_delete                  |               3        |                      3        |
| bulk_delete_normal           |               0.428571 |                      0.75     |
| damerau                      |               3        |                      4        |
| damerau_normal               |               0.428571 |                      1        |
| entropy                      |               0.101341 |                      0        |
| entropy_normal               |               0.101341 |                      0        |
| fast_birnbaum                |               0.666667 |                      0.7      |
| fast_birnbaum_normal         |               0.666667 |                      0.7      |
| fast_birnbaum_simil          |               7        |                      3        |
| fast_birnbaum_simil_normal   |               0.25     |                      0.3      |
| fragile_ends_simil           |               3        |                      3.5      |
| fragile_ends_simil_normal    |               0.5      |                      1        |
| jaccard                      |               0.7      |                      0        |
| jaccard_normal               |               0.7      |                      0        |
| jaro                         |               0.253968 |                      0.5      |
| jaro_normal                  |               0.253968 |                      0.5      |
| jaro_winkler                 |               0.253968 |                      0.5      |
| jaro_winkler_normal          |               0.253968 |                      0.5      |
| levenshtein                  |               3        |                      4        |
| levenshtein_normal           |               0.428571 |                      1        |
| mmcwpa                       |               0.538462 |                      0.387628 |
| mmcwpa_normal                |               0.538462 |                      0.387628 |
| ratcliff_obershelp           |               0.384615 |                      0.5      |
| ratcliff_obershelp_normal    |               0.384615 |                      0.5      |
| sorensen                     |               0.384615 |                      0        |
| sorensen_normal              |               0.384615 |                      0        |
| stemmatological_simil        |               3        |                      3        |
| stemmatological_simil_normal |               0.428571 |                      0.75     |
| subseq_jaccard               |               0.751556 |                      0.547008 |
| subseq_jaccard_normal        |               0.751556 |                      0.547008 |


## Changelog

Version 0.3:

  - Improvements to code quality, documentation, and references
  - Added new methods and scaffolding for future expansions

Version 0.2:

  - First release for new roadmap supporting sequences of any hashable Python
    datatype, importing code from other projects (mostly from `titivillus`)
    
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
> types. Version 0.3. Uppsala: Uppsala universitet.
> Available at: https://github.com/evotext/seqsim

In BibTeX:

```
@misc{Tresoldi2021seqsim,
  author = {Tresoldi, Tiago; Maurits, Luke; Dunn, Michael},
  title = {seqsim, a library for computing measures of distance and similarity for sequences of hashable data types. Version 0.3},
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
