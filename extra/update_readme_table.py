"""
Command-line script for updating the table of measures in the README.

The table is generated from the declarations of the measures, between the
`measures:start` and `measures:end` markers of `README.md`. Run it after
adding a measure or changing a declaration; `tests/test_docs.py` fails while
the README is out of date.
"""

import pathlib
import re

from seqsim._measure import markdown_table

README = pathlib.Path(__file__).parent.parent / "README.md"
MARKERS = re.compile(r"(<!-- measures:start.*?-->\n).*?(<!-- measures:end -->)", re.S)


def updated(text: str) -> str:
    """
    Returns the text of the README with the table of measures regenerated.
    """

    if not MARKERS.search(text):
        raise ValueError("The README has no measures:start/end markers.")
    return MARKERS.sub(lambda match: match[1] + markdown_table() + match[2], text)


def main():
    text = README.read_text(encoding="utf-8")
    README.write_text(updated(text), encoding="utf-8")


if __name__ == "__main__":
    main()
