"""
test_docs
=========

Runs the examples in the README and in the docstrings as doctests, so that
the documentation cannot drift from the code.
"""

# Import Python standard libraries
import doctest
import pathlib
import re

import pytest

# Import the library being tested
import seqsim
from seqsim import (
    _measure,
    alignment,
    common,
    compression,
    edit,
    ngrams,
    order,
    sequence,
    token,
)
from seqsim.tradition import _characters, _coverage, _export, _frame

ROOT = pathlib.Path(__file__).parent.parent
README = ROOT / "README.md"
DOCS = sorted(
    path
    for path in (ROOT / "docs").rglob("*.md")
    if path.name != "candidate_methods.md"
)


def _run(text, name, globs=None):
    parser = doctest.DocTestParser()
    globs = {"seqsim": seqsim} if globs is None else globs
    test = parser.get_doctest(text, globs, name, name, 0)
    runner = doctest.DocTestRunner(
        optionflags=doctest.ELLIPSIS | doctest.NORMALIZE_WHITESPACE
    )
    runner.run(test, clear_globs=False)

    # `DocTest` works on a copy of the namespace; propagate the names
    # defined in this block to the following ones
    globs.update(test.globs)
    results = runner.summarize(verbose=False)
    assert results.failed == 0, f"{results.failed} failed examples in {name}"
    return results.attempted


def test_readme_examples():
    blocks = re.findall(
        r"```python\n(.*?)```", README.read_text(encoding="utf-8"), re.S
    )
    assert sum(_run(block, "README") for block in blocks) > 0


@pytest.mark.parametrize("path", DOCS, ids=lambda path: path.name)
def test_documentation_examples(path, monkeypatch):
    # All blocks of a page share their namespace, as in a tutorial; examples
    # run from the folder with the data files, which readers download
    monkeypatch.chdir(ROOT / "docs" / "data")
    blocks = re.findall(r"```python\n(.*?)```", path.read_text(encoding="utf-8"), re.S)
    globs = {"seqsim": seqsim}
    for idx, block in enumerate(blocks):
        _run(block, f"{path.name}[{idx}]", globs)


@pytest.mark.parametrize(
    "module",
    [
        seqsim,
        _measure,
        alignment,
        common,
        compression,
        edit,
        ngrams,
        order,
        sequence,
        token,
        _frame,
        _coverage,
        _characters,
        _export,
    ],
)
def test_docstring_examples(module):
    finder = doctest.DocTestFinder()
    runner = doctest.DocTestRunner(optionflags=doctest.ELLIPSIS)
    for test in finder.find(module, globs={"seqsim": seqsim}):
        runner.run(test)
    assert runner.summarize(verbose=False).failed == 0
