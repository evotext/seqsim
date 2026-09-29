# Configuration file for the Sphinx documentation builder.
#
# See https://www.sphinx-doc.org/en/master/usage/configuration.html

import os
import sys

sys.path.insert(0, os.path.abspath(os.path.join("..", "src")))

import seqsim  # noqa: E402

# -- Project information -----------------------------------------------------

project = "seqsim"
copyright = "2021-2026, Tiago Tresoldi, Luke Maurits, Michael Dunn"
author = "Tiago Tresoldi, Luke Maurits, Michael Dunn"
release = seqsim.__version__
version = release

# -- General configuration ---------------------------------------------------

extensions = [
    "myst_parser",
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.mathjax",
    "sphinx.ext.viewcode",
    "sphinx_copybutton",
]

source_suffix = {".md": "markdown", ".rst": "restructuredtext"}
root_doc = "index"

myst_enable_extensions = [
    "colon_fence",
    "deflist",
    "dollarmath",
    "amsmath",
    "attrs_block",
]
myst_heading_anchors = 3

autodoc_member_order = "bysource"
autodoc_typehints = "description"

exclude_patterns = [
    "_build",
    "Thumbs.db",
    ".DS_Store",
    "candidate_methods.md",
    "data/README.md",
    "methods/_measures.md",
]

# -- Generated content -------------------------------------------------------

# The section of the method pages describing each measure of `distance()`;
# the table of measures is generated from their declarations, so that it
# cannot drift from the code
SECTIONS = {
    "levenshtein": "levenshtein",
    "damerau": "damerau-levenshtein",
    "osa": "optimal-string-alignment",
    "indel": "indel-and-lcs",
    "lcs": "indel-and-lcs",
    "levenshtein_gld": "normalized-edit-distances",
    "damerau_gld": "normalized-edit-distances",
    "indel_gld": "normalized-edit-distances",
    "levenshtein_ned": "normalized-edit-distances",
    "bulk_delete": "bulk-delete",
    "fragile_ends": "fragile-ends",
    "stemmatological": "stemmatological",
    "block_move": "block-moves",
    "gst": "greedy-string-tiling",
    "jaro": "jaro-and-jaro-winkler",
    "jaro_winkler": "jaro-and-jaro-winkler",
    "mmcwpa": "mmcwpa",
    "birnbaum": "birnbaum",
    "ulam": "ulam",
    "kendall_tau": "kendall-tau",
    "footrule": "spearman-footrule",
    "cayley": "cayley",
    "block_interchange": "block-interchange",
    "breakpoint": "breakpoints",
    "nw": "global-alignment",
    "jaccard": "jaccard",
    "sorensen": "sørensen-dice",
    "subseq_jaccard": "sub-sequence-jaccard",
    "qgram": "q-grams",
    "ratcliff_obershelp": "ratcliff-obershelp",
    "entropy_ncd": "entropy-ncd",
    "lzma_ncd": "lzma-ncd",
    "lz76": "lempel-ziv",
}


def _section(info):
    module = info.name.split(".")[0]
    return f"{module}.md#{SECTIONS[info.key]}"


_table = seqsim._measure.markdown_table(link=_section)
_target = os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "methods", "_measures.md"
)
with open(_target, "w", encoding="utf-8") as handler:
    handler.write(_table)

# -- Options for HTML output -------------------------------------------------

html_theme = "furo"
html_title = f"seqsim {release}"
html_static_path = ["_static"]
html_logo = "scriptorium_small.jpg"
html_theme_options = {
    "source_repository": "https://github.com/evotext/seqsim",
    "source_branch": "main",
    "source_directory": "docs/",
}

# Copy only the code in examples, not the prompts or outputs
copybutton_prompt_text = r">>> |\.\.\. "
copybutton_prompt_is_regexp = True
