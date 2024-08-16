import sys
from pathlib import Path

import plinder

# TODO: Use proper __version__ attribute from the Plinder package
plinder.__version__ = "0.1.0"

DOC_PATH = Path(__file__).parent
PACKAGE_PATH = DOC_PATH.parent / "src"
COLUMN_REFERENCE_PATH = DOC_PATH.parent / "column_descriptions"

# Include documentation in PYTHONPATH
# in order to import modules for API doc generation etc.
sys.path.insert(0, str(DOC_PATH))
import tablegen

# Pregeneration of files
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
    # "sphinx.ext.linkcode",
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

templates_path = ["templates"]
source_suffix = {
    '.rst': 'restructuredtext',
    '.ipynb': 'myst-nb',
    '.myst': 'myst-nb',
}
master_doc = "index"

project = "PLINDER"
copyright = "2024, PLINDER Development Team"
author = "Plinder Development Team"
version = plinder.__version__
release = plinder.__version__

pygments_style = "sphinx"

#### HTML ####

html_theme = "pydata_sphinx_theme"

html_static_path = ["static"]
html_css_files = [
    "plinder.css",
    # For rendering the column descriptions table
    "https://cdn.datatables.net/2.1.3/css/dataTables.dataTables.css",
    # Get fonts from Google Fonts CDN
    "https://fonts.googleapis.com/css2"
    "?family=Geologica:wght@100..900"
    "&family=Montserrat:ital,wght@0,100..900;1,100..900"
    "&display=swap",
]
html_js_files = [
    # For rendering the column descriptions table
    "https://code.jquery.com/jquery-3.7.1.js",
    "https://cdn.datatables.net/2.1.3/js/dataTables.js",
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
}
html_context = {
    "github_user": "plinder-org",
    "github_repo": "plinder",
    "github_version": "main",
    "doc_path": "doc",
}


#### App setup ####


# def setup(app):
#    app.connect("autodoc-skip-member", apidoc.skip_nonrelevant)
