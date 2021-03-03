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

The library offers different methods to compare sequences of hashable Python elements
(not only characters in a string).

```python
>> import seqsim
>> test1_seq_a = "kitten"
>> test1_seq_b = "sitting"
>> test2_seq_a = (1, 2, 3)
>> test2_seq_b = [1, 2, 3]
>> test3_seq_a = (1, 2, 3, 4, 5)
>> test3_seq_b = (1, 2, 4, 3, 6, 7)
>> test4_seq_a = (1, 2, 3)
>> test4_seq_b = ["a", "b", "c", "d"]
>> seqsim.edit.prev_edit_distance(test1_seq_a, test1_seq_b)
3.0
>> seqsim.edit.jaccard_dist(test3_seq_a, test3_seq_b)
0.4285714285714286
>> seqsim.edit.mmcwpa_dist(test3_seq_a, test3_seq_b)
0.5546382285848768
```

Full documentation will be offered in future releases. For the moment, the library
usage is illustrated by set of tests in the `tests/` directory.


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

> Tresoldi, Tiago; Maurits, Luke; Dunn, Michael. (2021). seqsim, a library for computing
> measures of similarity for sequences of hashable data types. Version 0.3.
> Uppsala: Uppsala universitet. Available at: https://github.com/evotext/seqsim

In BibTeX:

```
@misc{Tresoldi2021seqsim,
  author = {Tresoldi, Tiago; Maurits, Luke; Dunn, Michael},
  title = {seqsim, a library for computing measures of similarity for sequences of hashable data types. Version 0.3},
  howpublished = {\url{https://github.com/evotext/seqsim}},
  address = {Uppsala},
  publisher = {Uppsala universitet},
  year = {2021},
}

## References

The image at the top of this file is derived from Yves de Saint-Denis, *Vie et martyre de saint
Denis et de ses compagnons, versions latine et française*. It is available in high
resolution from [Bibliothèque nationale de France, Département des Manuscrits, Français 2090,
fol. 12v.](http://gallica.bnf.fr/ark:/12148/btv1b8447296x/f30.item)

References to the various implementation are available in the source code comments and in
the [online documentation](https://seqsim.readthedocs.io/en/latest/?badge=latest).