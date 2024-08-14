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
import apidoc
import tablegen
import viewcode

# Pregeneration of files
#apidoc.create_api_doc(PACKAGE_PATH, DOC_PATH / "apidoc")
#tablegen.generate_table(COLUMN_REFERENCE_PATH, DOC_PATH / "table.md")

#### Source code link ###

#linkcode_resolve = viewcode.linkcode_resolve

#### General ####

extensions = [
    "jupyter_sphinx",
    "sphinx.ext.autodoc",
    #"sphinx.ext.autosummary",
    "sphinx.ext.doctest",
    "sphinx.ext.mathjax",
    #"sphinx.ext.linkcode",
    "sphinx.ext.todo",
    "sphinx_design",
    "sphinx_copybutton",
    "numpydoc",
    "myst_parser",
    #"myst_nb",
]

# TODO: remove this for production
todo_include_todos = True

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
source_suffix = {".md": "markdown"}
master_doc = "index"

project = "Plinder"
copyright = '2024, Plinder Development Team'
author = 'Plinder Development Team'
version = plinder.__version__
release = plinder.__version__

pygments_style = "sphinx"

#### HTML ####

html_theme = "pydata_sphinx_theme"

html_static_path = ["static"]
html_css_files = ["plinder.css", "fonts.css"]
html_title = "Plinder"
html_logo = "static/assets/general/plinder_logo.webp"
html_favicon = "static/assets/general/plinder_icon.png"
html_baseurl = f"https://plinder-org.github.io/plinder/"
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
            "name": "PyPI",
            "url": "https://pypi.org/project/plinder/",
            "icon": "fa-solid fa-box-open",
            "type": "fontawesome",
        },
    ],
    "use_edit_page_button": True,
    "show_prev_next": False,
    "show_toc_level": 2,
}
html_sidebars = {
    # No primary sidebar for these pages
    "tutorial/dataset": [],
    "tutorial/api": [],
}
html_context = {
    "github_user": "plinder-org",
    "github_repo": "plinder",
    "github_version": "main",
    "doc_path": "doc",
}


#### App setup ####


#def setup(app):
#    app.connect("autodoc-skip-member", apidoc.skip_nonrelevant)
