# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'UKB PPP Download'
copyright = '2026, Natacha Galmiche'
author = 'Natacha Galmiche'

# The full version, including alpha/beta/rc tags
release = '0.1.1'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

# Note that as of December 2023, if 'm2r' and 'myst_parser'
# can not be activated simultaneously, however m2r seems to include
# myst_parser
extensions = [
    # To generate documentation from docstrings
    'sphinx.ext.autodoc',
    # To be able to use numpy styles or google styles of docstrings
    'sphinx.ext.napoleon',
    # To get only the first line of all functions and potentially use
    # it as a toctree
    'sphinx.ext.autosummary',
    # To tell sphinx that we are also using markdown
    # 'myst_parser',
    # To be able to include md files directly in rst files
    # adds the mdinclude directive
    'm2r',
]

# Uncomment if m2r is NOT included in the extension
# But myst_parser IS
# To tell sphinx to look also for .md files.
# source_suffix = ['.rst', '.md']

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store', "../scripts", "../data"]



# -- Options for HTML output -------------------------------------------------
# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

# -- Autodoc options -------------------------------------------------

# We don't want to document tests.
autodoc_mock_imports = [
    "ukbppp_dl.tests", "scripts",
]

autodoc_default_options = {
    # will document all class member methods and properties.
    'members': True,
    'member-order': 'bysource',
    # will also generate document for the special member (__foo__)
    'special-members': '__call__',
    # will also generate document for the members not having docstrings:
    'undoc-members': True,
    # will also generate document for private members (_ or __)
    'private_members': False,
    'exclude-members': '__weakref__'
}
