# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'plinder'
copyright = '2024, Plinder Development Team'
author = 'Plinder Development Team'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.doctest",
    "sphinx.ext.autosummary",
    "sphinx.ext.todo",
    "sphinx.ext.coverage",
    "sphinx.ext.mathjax",
    "sphinx.ext.viewcode",
    "sphinx.ext.napoleon",
    "sphinx.ext.githubpages",
    "myst_nb",
]

source_suffix = {
    ".rst": "restructuredtext",
    ".ipynb": "myst-nb",
    ".myst": "myst-nb",
    ".md": "myst-nb",
}

nb_custom_formats = {".ipynb": ["jupytext.reads", {"fmt": "ipynb"}]}
nb_execution_timeout = 720
nb_kernel_rgx_aliases = {"plinder.*": "python3"}

templates_path = ['_templates']
exclude_patterns = ["build", ".DS_Store"]

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


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'piccolo_theme'
html_static_path = ['_static']
html_logo = "../assets/plinder.png"
html_theme_options = {
    "source_url": "https://github.com/plinder-org/plinder/",
    "source_icon": "github",
    "show_theme_credit": False,
}
html_favicon = "_static/favicon.png"
htmlhelp_basename = "PlinderDoc"
