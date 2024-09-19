import os
import sys
from pathlib import Path

import plinder

DOC_PATH = Path(__file__).parent
COLUMN_REFERENCE_PATH = DOC_PATH.parent / "src" / "plinder" / "data" / "column_descriptions"

# Avoid verbose logs in rendered notebooks
os.environ["PLINDER_LOG_LEVEL"] = "0"

# Include documentation in PYTHONPATH
# in order to import modules for API doc generation etc.
sys.path.insert(0, str(DOC_PATH))
import apidoc
import tablegen
import viewcode

# Pregeneration of files
for package in ["plinder.core", "plinder.core.scores", "plinder.core.loader"]:
    apidoc.generate_api_reference(package, DOC_PATH / "api" / package.split(".")[-1])
tablegen.generate_table(COLUMN_REFERENCE_PATH, DOC_PATH / "table.html")

#### Source code link ###

# linkcode_resolve = viewcode.linkcode_resolve

#### General ####

extensions = [
    "jupyter_sphinx",
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.doctest",
    "sphinx.ext.mathjax",
    "sphinx.ext.linkcode",
    "sphinx.ext.todo",
    "sphinx_design",
    "sphinx_copybutton",
    "numpydoc",
    "myst_nb",
]

nb_custom_formats = {".ipynb": ["jupytext.reads", {"fmt": "ipynb"}]}
nb_execution_timeout = 720
nb_kernel_rgx_aliases = {"plinder.*": "python3"}
myst_enable_extensions = [
    "amsmath",
    "attrs_inline",
    "colon_fence",
    "deflist",
    "dollarmath",
    "fieldlist",
    "html_admonition",
    "html_image",
    "linkify",
    "replacements",
    "smartquotes",
    "strikethrough",
    "substitution",
    "tasklist",
]
myst_url_schemes = ("http", "https", "mailto")

numpydoc_show_class_members = False
# Prevent autosummary from using sphinx-autogen, since it would
# overwrite the document structure given by apidoc.json
autosummary_generate = False
linkcode_resolve = viewcode.linkcode_resolve

templates_path = ["templates"]
source_suffix = {
    ".rst": "restructuredtext",
    ".ipynb": "myst-nb",
    ".myst": "myst-nb",
}
master_doc = "index"

project = "PLINDER"
copyright = "2024, PLINDER Development Team"
author = "PLINDER Development Team"
version = plinder.__version__
release = plinder.__version__

pygments_style = "sphinx"

#### HTML ####

html_theme = "pydata_sphinx_theme"

html_static_path = ["static"]
html_css_files = [
    "plinder.css",
    # Get fonts from Google Fonts CDN
    "https://fonts.googleapis.com/css2"
    "?family=Geologica:wght@100..900"
    "&family=Montserrat:ital,wght@0,100..900;1,100..900"
    "&display=swap",
]
html_title = "PLINDER"
html_logo = "static/assets/general/plinder_logo.png"
html_favicon = "static/assets/general/plinder_icon.png"
html_baseurl = "https://plinder-org.github.io/plinder/"
html_theme_options = {
    "header_links_before_dropdown": 7,
    "pygments_light_style": "friendly",
    "icon_links": [
        {
            "name": "GitHub",
            "url": "https://github.com/plinder-org/plinder/",
            "icon": "fa-brands fa-github",
            "type": "fontawesome",
        },
        {
            "name": "Article",
            "url": "https://www.biorxiv.org/content/10.1101/2024.07.17.603955v3",
            "icon": "fa-solid fa-file-lines",
            "type": "fontawesome",
        },
    ],
    "use_edit_page_button": True,
    "show_prev_next": False,
    "show_toc_level": 2,
}
html_sidebars = {
    # No primary sidebar for these pages
    "dataset": [],
    "evaluation": [],
    "citation": [],
}
html_context = {
    "github_user": "plinder-org",
    "github_repo": "plinder",
    "github_version": "main",
    "doc_path": "doc",
}
html_scaled_image_link = False


#### App setup ####


def setup(app):
    app.connect("autodoc-skip-member", apidoc.skip_nonrelevant)
