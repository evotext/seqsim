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

exclude_patterns = ["_build", "Thumbs.db", ".DS_Store", "candidate_methods.md", "data/README.md"]

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
