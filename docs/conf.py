#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# RADIS documentation build configuration file, created by
# sphinx-quickstart on Tue Feb  6 03:32:15 2018.
#
# This file is execfile()d with the current directory set to its
# containing dir.
#
# Note that not all possible configuration values are present in this
# autogenerated file.
#
# All configuration values have a default; values that are commented out
# serve to show the default.
#
# --------
# For help on Sphinx and readthedocs see:
# https://samnicholls.net/2016/06/15/how-to-sphinx-readthedocs/
# --------

import os
import sys

import sphinx_gallery.gen_rst

# %% Custom Example header
# https://github.com/sphinx-gallery/sphinx-gallery/issues/978
sphinx_gallery.gen_rst.EXAMPLE_HEADER = """
.. DO NOT EDIT.
.. THIS FILE WAS AUTOMATICALLY GENERATED BY SPHINX-GALLERY.
.. TO MAKE CHANGES, EDIT THE SOURCE PYTHON FILE:
.. "{0}"
.. LINE NUMBERS ARE GIVEN BELOW.

.. only:: html

    .. note::
        :class: sphx-glr-download-link-note

        Run this example online :

        - Click :ref:`here <sphx_glr_download_{1}>`
          to download the full example code{2}

        - Then start `Radis-Lab <https://radis.github.io/radis-lab/>`__,
          upload the Jupyter notebook, and run it from there.

.. rst-class:: sphx-glr-example-title

.. _sphx_glr_{1}:

"""

#%%
# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
sys.path.insert(0, os.path.abspath(".."))


# %% ------------------------------------
# Added EP 2018:
# Auto-generate files with sphinx.apidoc
# (else it requires /docs/source files to be generated manually and committed
# to the git directory)
#
# Reference:
# https://github.com/rtfd/readthedocs.org/issues/1139
#

# %% --------------------------------------


# -- General configuration ------------------------------------------------

# If your documentation needs a minimal Sphinx version, state it here.
# needs_sphinx = '1.0'

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.doctest",
    "sphinx.ext.mathjax",
    #'sphinx.ext.viewcode',
    "sphinx.ext.autosummary",
    "sphinx_gallery.gen_gallery",
    #'numpydoc',
    #'sphinxcontrib.napoleon',
    "sphinx.ext.napoleon",
    "sphinx_autodoc_defaultargs",
    "sphinx.ext.intersphinx",
    "sphinx.ext.inheritance_diagram",
    "sphinxcontrib.apidoc",
    "sphinx.ext.linkcode",
]

sphinx_gallery_conf = {
    "examples_dirs": "../examples",  # path to your example scripts
    "gallery_dirs": "auto_examples",  # path to where to save gallery generated output
    # to make references clickable
    "doc_module": "radis",
    "reference_url": {
        "radis": None,
    },
    # directory where function/class granular galleries are stored
    "backreferences_dir": "source/backreferences",
    # Modules for which function/class level galleries are created.
    "doc_module": ("radis"),
    "inspect_global_variables": True,
    "show_signature": False,
}


def linkcode_resolve(domain, info):
    """for sphinx.ext.linkcode"""
    if domain != "py":
        return None
    if not info["module"]:
        return None
    filename = info["module"].replace(".", "/")
    return "https://github.com/radis/radis/tree/develop/%s.py" % filename


# used to mini-galleries : https://sphinx-gallery.github.io/stable/configuration.html#add-mini-galleries-for-api-documentation
autosummary_generate = True

# %% ------------------------------------
# Added EP 2018:
# Auto-generate files with sphinx.apidoc
# (else it requires /docs/source files to be generated manually and committed
# to the git directory)
#
# Reference:
# https://github.com/rtfd/readthedocs.org/issues/1139
#


def run_apidoc(_):

    argv = [
        "-f",
        "-e",
        "-o",
        "source",
        "--separate",
        "../radis",
    ]

    try:
        # Sphinx 1.7+
        from sphinx.ext import apidoc

        apidoc.main(argv)
    except ImportError:
        # Sphinx 1.6 (and earlier)
        from sphinx import apidoc

        argv.insert(0, apidoc.__file__)
        apidoc.main(argv)


def setup(app):
    app.connect("builder-inited", run_apidoc)
    app.add_css_file("custom.css")  #  for scrollable sidebar


# %%

# Reference other packages
intersphinx_mapping = {
    "astropy": ("https://docs.astropy.org/en/stable/", None),
    "astroquery": ("https://astroquery.readthedocs.io/en/latest/", None),
    "cantera": ("https://www.cantera.org/documentation/docs-2.6/sphinx/", None),
    "fitroom": ("https://fitroom.readthedocs.io/en/latest/", None),
    "habanero": ("https://habanero.readthedocs.io/en/latest/", None),
    "joblib": ("https://joblib.readthedocs.io/en/latest/", None),
    "lmfit": ("https://lmfit.github.io/lmfit-py/", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
    "matplotlib": ("https://matplotlib.org/stable/", None),
    "pandas": ("https://pandas.pydata.org/pandas-docs/stable/", None),
    "pytexit": ("https://pytexit.readthedocs.io/en/latest/", None),
    "scipy": ("https://docs.scipy.org/doc/scipy/", None),
    "seaborn": ("https://seaborn.pydata.org/", None),
    "specutils": ("https://specutils.readthedocs.io/en/stable/", None),
    "vaex": ("https://vaex.readthedocs.io/en/latest/", None),
}

napoleon_google_docstring = False
napoleon_use_param = False  # fixes https://github.com/sphinx-doc/sphinx/issues/2227

# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]

# The suffix(es) of source filenames.
# You can specify multiple suffix as a list of string:
# source_suffix = ['.rst', '.md']
source_suffix = ".rst"

# The encoding of source files.
# source_encoding = 'utf-8-sig'

# The master toctree document.
master_doc = "index"

# General information about the project.
project = "RADIS"
copyright = "2021, Erwan Pannier and the 🌱 RADIS contributors (https://github.com/radis/radis/graphs/contributors)"
author = (
    "Erwan Pannier, Dirk van den Bekerom, Nicolas Minesi, "
    + "et al. (https://github.com/radis/radis/graphs/contributors)"
)

# The version info for the project you're documenting, acts as a replacement for
# |version| and |release|, also used in various other places throughout the
# built documents.
#

with open("../radis/__version__.txt") as version_file:
    __version__ = version_file.read().strip()

# The short X.Y version.
version = __version__
# The full version, including alpha/beta/rc tags.
release = __version__

# The language for content autogenerated by Sphinx. Refer to documentation
# for a list of supported languages.
#
# This is also used if you do content translation via gettext catalogs.
# Usually you set "language" from the command line for these cases.
language = "en"

# There are two options for replacing |today|: either, you set today to some
# non-false value, then it is used:
# today = ''
# Else, today_fmt is used as the format for a strftime call.
# today_fmt = '%B %d, %Y'

# List of patterns, relative to the source directory, that match files and
# directories to ignore when looking for source files.
exclude_patterns = ["_build", "dev/_*", "lbl/_*", "spectrum/_*", "references/_*"]

# The reST default role (used for this markup: `text`) to use for all
# documents.
# default_role = None
default_role = "py:obj"  # fixes https://github.com/spacetelescope/specviz/pull/77

# If true, '()' will be appended to :func: etc. cross-reference text.
# add_function_parentheses = True

# If true, the current module name will be prepended to all description
# unit titles (such as .. function::).
add_module_names = False

# If true, sectionauthor and moduleauthor directives will be shown in the
# output. They are ignored by default.
# show_authors = False

# The name of the Pygments (syntax highlighting) style to use.
pygments_style = "sphinx"

# A list of ignored prefixes for module index sorting.
# modindex_common_prefix = []

# If true, keep warnings as "system message" paragraphs in the built documents.
# keep_warnings = False

# If true, `todo` and `todoList` produce output, else they produce nothing.
todo_include_todos = False


# -- Options for HTML output ----------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
html_theme = "alabaster"  #'bizstyle' #'nature'

# Theme options are theme-specific and customize the look and feel of a theme
# further.  For a list of options available for each theme, see the
# documentation.
html_theme_options = {
    "description": "Radiative Solver",
    #'logo': 'radis_ico.png',
    "logo_name": True,
    "github_user": "radis",
    "github_repo": "radis",
    "github_button": True,
    "github_type": "star",
    "github_banner": False,
    "travis_button": False,
    "codecov_button": False,
    "sidebar_includehidden": True,
    "fixed_sidebar": True,
    "analytics_id": "UA-113616205-1",
    "link": "#7A306C",
    "extra_nav_links": {
        "RADIS Website": "https://radis.github.io/",
        "Video Tutorials": "https://www.youtube.com/channel/UCO-7NXkubTAiGGxXmvtQlsA",
    },
}

# Add any paths that contain custom themes here, relative to this directory.
# html_theme_path = []

# The name for this set of Sphinx documents.  If None, it defaults to
# "<project> v<release> documentation".
# html_title = None

# A shorter title for the navigation bar.  Default is the same as html_title.
# html_short_title = None

# The name of an image file (relative to this directory) to place at the top
# of the sidebar.
html_logo = "radis_ico.png"

# The name of an image file (within the static path) to use as the favicon of the
# docs.  This file should be a Windows icon file (.ico) being 16x16 or 32x32
# pixels large.
# html_favicon = None

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the built-in static files,
# so a file named "default.css" will overwrite the built-in "default.css".
html_static_path = ["_static"]

# Add any extra paths that contain custom files (such as robots.txt or
# .htaccess) here, relative to this directory. These files are copied
# directly to the root of the documentation.
# html_extra_path = []

# If not '', a 'Last updated on:' timestamp is inserted at every page bottom,
# using the given strftime format.
# html_last_updated_fmt = '%b %d, %Y'

# If true, SmartyPants will be used to convert quotes and dashes to
# typographically correct entities.
# html_use_smartypants = True

# Custom sidebar templates, maps document names to template names.
html_sidebars = {
    "**": [
        "about.html",
        "navigation.html",
        #'relations.html',
        #'description.html',
        #'example.html',
        "searchbox.html",
        "buttons.html",
    ]
}


# The default values of all documented arguments, and undocumented arguments if enabled, are automatically detected and added to the docstring.
# It also detects existing documentation of default arguments with the text unchanged.

rst_prolog = (
    """
.. |default| raw:: html

    <div class="default-value-section">"""
    + ' <span class="default-value-label">Default:</span>'
)
# Additional templates that should be rendered to pages, maps page names to
# template names.
# html_additional_pages = {}

# If false, no module index is generated.
# html_domain_indices = True

# If false, no index is generated.
# html_use_index = True

# If true, the index is split into individual pages for each letter.
# html_split_index = False

# If true, links to the reST sources are added to the pages.
# html_show_sourcelink = True

# If true, "Created using Sphinx" is shown in the HTML footer. Default is True.
# html_show_sphinx = True

# If true, "(C) Copyright ..." is shown in the HTML footer. Default is True.
# html_show_copyright = True

# If true, an OpenSearch description file will be output, and all pages will
# contain a <link> tag referring to it.  The value of this option must be the
# base URL from which the finished HTML is served.
# html_use_opensearch = ''

# This is the file name suffix for HTML files (e.g. ".xhtml").
# html_file_suffix = None

# Language to be used for generating the HTML full-text search index.
# Sphinx supports the following languages:
#   'da', 'de', 'en', 'es', 'fi', 'fr', 'h', 'it', 'ja'
#   'nl', 'no', 'pt', 'ro', 'r', 'sv', 'tr'
# html_search_language = 'en'

# A dictionary with options for the search language support, empty by default.
# Now only 'ja' uses this config value
# html_search_options = {'type': 'default'}

# The name of a javascript file (relative to the configuration directory) that
# implements a search results scorer. If empty, the default will be used.
# html_search_scorer = 'scorer.js'

# Output file base name for HTML help builder.
htmlhelp_basename = "RADISdoc"

# -- Options for LaTeX output ---------------------------------------------

latex_elements = {
    # The paper size ('letterpaper' or 'a4paper').
    #'papersize': 'letterpaper',
    # The font size ('10pt', '11pt' or '12pt').
    #'pointsize': '10pt',
    # Additional stuff for the LaTeX preamble.
    "preamble": "\setcounter{tocdepth}{3}",
    # Latex figure (float) alignment
    #'figure_align': 'htbp',
}

# Grouping the document tree into LaTeX files. List of tuples
# (source start file, target name, title,
#  author, documentclass [howto, manual, or own class]).
latex_documents = [
    (master_doc, "RADIS.tex", "RADIS Documentation", author, "manual"),
]

# The name of an image file (relative to this directory) to place at the top of
# the title page.
# latex_logo = None

# For "manual" documents, if this is true, then toplevel headings are parts,
# not chapters.
latex_use_parts = True

# If true, show page references after internal links.
# latex_show_pagerefs = False

# If true, show URL addresses after external links.
# latex_show_urls = False

# Documents to append as an appendix to all manuals.
# latex_appendices = []

# If false, no module index is generated.
# latex_domain_indices = True


# -- Options for manual page output ---------------------------------------

# One entry per manual page. List of tuples
# (source start file, name, description, authors, manual section).
man_pages = [(master_doc, "radis", "RADIS Documentation", [author], 1)]

# If true, show URL addresses after external links.
# man_show_urls = False


# -- Options for Texinfo output -------------------------------------------

# Grouping the document tree into Texinfo files. List of tuples
# (source start file, target name, title, author,
#  dir menu entry, description, category)
texinfo_documents = [
    (
        master_doc,
        "RADIS",
        "RADIS Documentation",
        author,
        "RADIS",
        "A fast line-by-line code for high resolution infrared molecular spectra: https://radis.github.io/",
        "Miscellaneous",
    ),
]

# Documents to append as an appendix to all manuals.
# texinfo_appendices = []

# If false, no module index is generated.
# texinfo_domain_indices = True

# How to display URL addresses: 'footnote', 'no', or 'inline'.
# texinfo_show_urls = 'footnote'

# If true, do not generate a @detailmenu in the "Top" node's menu.
# texinfo_no_detailmenu = False


# Extra
numpydoc_show_class_members = False
