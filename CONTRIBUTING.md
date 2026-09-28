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

- Methods take two sequences of arbitrary hashable elements as positional
  arguments; all other parameters, including `normal`, are keyword-only.
- Name the function after its properties: `_dist` only for true metrics,
  `_dissim` for other dissimilarities (`0.0` for identical sequences), and
  `_simil` for similarities.
- Return a float, handle empty sequences (two empty sequences score `0.0`),
  and make sure `normal=True` returns values in range [0..1].
- Add the method to `METHODS` in `src/seqsim/__init__.py` and to the lists in
  `tests/test_properties.py`, which check the promised properties. Expected
  values in the tests should be verified by hand or against an independent
  reference.
- Regenerate the README table with `python extra/readme_compare.py`, and
  document the method in the corresponding page of `docs/methods`.

## Code of conduct

Be respectful and constructive. Harassment or discriminatory behaviour of any
kind is not tolerated.
