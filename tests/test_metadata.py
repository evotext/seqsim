"""Consistency of the version stated in the different project files."""

import re
from pathlib import Path

import seqsim

ROOT = Path(__file__).resolve().parent.parent


def test_citation_version():
    text = (ROOT / "CITATION.cff").read_text(encoding="utf-8")
    match = re.search(r"^version:\s*(\S+)", text, re.M)
    assert match and match.group(1) == seqsim.__version__


def test_changelog_mentions_version():
    text = (ROOT / "CHANGELOG.md").read_text(encoding="utf-8")
    assert f"## Version {seqsim.__version__}" in text


def test_typed_marker():
    assert (Path(seqsim.__file__).parent / "py.typed").is_file()
