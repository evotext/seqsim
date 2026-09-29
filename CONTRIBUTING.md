# Contributing to seqsim

Contributions are welcome, from bug reports to new methods. Please use the
GitHub [issues](https://github.com/evotext/seqsim/issues) and pull requests.

## Development setup

```bash
git clone https://github.com/evotext/seqsim.git
cd seqsim
python -m pip install -e ".[dev]"
```

Before submitting a pull request, please run:

```bash
black src tests
flake8 src tests --max-line-length=127 --extend-ignore=E203
pytest
```

## Documentation

The documentation is written in Markdown (MyST) in the `docs` folder. To
build it locally:

```bash
python -m pip install -e ".[docs]"
cd docs
sphinx-build -b html . _build/html
```

All examples in the documentation (lines starting with `>>>`) are run by the
test suite, so they must produce exactly the output shown.

## Adding a method

A method is a function decorated with `seqsim._measure.measure`, which
declares its kind and its claims and applies the conventions shared by all
measures (see `CONTEXT.md` for the vocabulary):

```python
@measure(key="levenshtein", kind="dist", triangle="yes", bound="max_len")
def levenshtein_dist(seq_x, seq_y):
    ...
```

- The body takes the two sequences (and the keyword-only options of the
  method, but not `normal`) and returns the raw value for that order of the
  arguments. The decorator validates the options (`check=`), applies the
  empty rule (`empty="max"`), symmetrizes (`symmetrize=`), normalizes
  (`bound=`), returns a float, and documents `normal`.
- Name the function after its properties: `_dist` only for true metrics,
  `_dissim` for other dissimilarities (`0.0` for identical sequences), and
  `_simil` for similarities.
- Declare the claims: `identity` and `triangle` are `"yes"`, `"no"` (with a
  counterexample), `"unproven"`, or `"conditional"` (with a condition).
  `tests/test_properties.py` tests every declared claim and verifies every
  counterexample, with no changes needed there. Expected values in other
  tests should be verified by hand or against an independent reference.
- With a `key`, the method is available in `distance()` and `METHODS`.
  Regenerate the table of the README with
  `python extra/update_readme_table.py` (the table of the documentation is
  generated when it is built, but the section of the method must be added to
  `SECTIONS` in `docs/conf.py`), and document the method in the
  corresponding page of `docs/methods`. The snapshot of results
  (`tests/data/snapshot.json`) is regenerated with
  `python tests/snapshot_cases.py`, only for intended changes.

## Code of conduct

Be respectful and constructive. Harassment or discriminatory behaviour of any
kind is not tolerated.
