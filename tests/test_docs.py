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
from seqsim import alignment, common, compression, edit, ngrams, order, sequence, token

README = pathlib.Path(__file__).parent.parent / "README.md"


def _run(text, name):
    parser = doctest.DocTestParser()
    test = parser.get_doctest(text, {"seqsim": seqsim}, name, str(README), 0)
    runner = doctest.DocTestRunner(optionflags=doctest.ELLIPSIS)
    runner.run(test)
    results = runner.summarize(verbose=False)
    assert results.failed == 0
    return results.attempted


def test_readme_examples():
    blocks = re.findall(
        r"```python\n(.*?)```", README.read_text(encoding="utf-8"), re.S
    )
    assert sum(_run(block, "README") for block in blocks) > 0


@pytest.mark.parametrize(
    "module",
    [seqsim, alignment, common, compression, edit, ngrams, order, sequence, token],
)
def test_docstring_examples(module):
    finder = doctest.DocTestFinder()
    runner = doctest.DocTestRunner(optionflags=doctest.ELLIPSIS)
    for test in finder.find(module, globs={"seqsim": seqsim}):
        runner.run(test)
    assert runner.summarize(verbose=False).failed == 0
