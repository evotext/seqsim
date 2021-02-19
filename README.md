# seqsim

[![PyPI](https://img.shields.io/pypi/v/seqsim.svg)](https://pypi.org/project/seqsim)
[![CI](https://github.com/tresoldi/seqsim/actions/workflows/main.yml/badge.svg)](https://github.com/tresoldi/seqsim/actions/workflows/main.yml)

Python library for computing measures of similarity for sequences of hashable data types

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
>> seqsim.distance.edit_distance(test1_seq_a, test1_seq_b)
3.0
>> seqsim.distance.jaccard_distance(test3_seq_a,test3_seq_b)
0.4285714285714286
>> seqsim.distance.mmcwpa_distance(test3_seq_a, test3_seq_b)
0.5546382285848768
```

Full documentation will be offered in future releases. For the moment, the library
usage is illustrated by set of tests in the `tests/` directory.


## Changelog

Version 0.2:

  - First release for new roadmap supporting sequences of any hashable Python
    datatype, importing code from other projects (mostly from `titivillus`)
    
## Community guidelines

While the author can be contacted directly for support, it is recommended that third 
parties use GitHub standard features, such as issues and pull requests, to contribute, 
report problems, or seek support.

Contributing guidelines, including a code of conduct, can be found in the
`CONTRIBUTING.md` file.

## Author and citation

The library is developed by Tiago Tresoldi (tiago.tresoldi@lingfil.uu.se).

If you use `seqsim`, please cite it as:

> Tresoldi, Tiago (2021). seqsim, a library for computing measures of similarity for
> sequences of hashable data types. Version 0.2. Uppsala: Uppsala Universitet.
> Available at: https://github.com/tresoldi/seqsim

In BibTeX:

```
@misc{Tresoldi2021titivillus,
  author = {Tresoldi, Tiago},
  title = {seqsim, a library for computing measures of similarity for sequences of hashable data types},
  howpublished = {\url{https://github.com/tresoldi/seqsim}},
  address = {Uppsala},
  publisher = {Uppsala Universitet},
  year = {2021},
}
