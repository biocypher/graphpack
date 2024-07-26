# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import os
import sys

sys.path.insert(0, os.path.abspath('.'))
sys.path.insert(0, os.path.abspath(os.path.join('..', '..')))
sys.path.insert(0, os.path.abspath('../../src/graphpack'))
sys.path.insert(0, os.path.abspath('../../src/graphpack/demo'))

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'graphpack'
copyright = '2024, Daniele'
author = 'Daniele Bottazzi'
version = '0.1.2'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.autosectionlabel',
    'sphinx.ext.napoleon',
    'sphinx.ext.autosummary',
    'sphinx.ext.doctest',
    'sphinx.ext.viewcode',
    'sphinx.ext.mathjax',
    'sphinx_togglebutton',
    'sphinxcontrib.bibtex',
]

autodoc_default_options = {
    'members': True,
    'undoc-members': True,
    'private-members': True,
    'special-members': '__init__',
    'inherited-members': True,
    'show-inheritance': True,
}

bibtex_bibfiles = ['refs.bib']
bibtex_default_style = 'plain'
templates_path = ['_templates']
exclude_patterns = []

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_theme_options = {
    "sidebarwidth": '25%',
}

templates_path = ['_templates']

html_static_path = ['_static']
html_extra_path = ['_files']

togglebutton_hint = "Click to see the algorithm pseudocode, parameters and hints"
togglebutton_hint_hide = "Click to hide"
