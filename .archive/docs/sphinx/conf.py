# -*- coding: utf-8 -*-
#
# Configuration file for the Sphinx documentation builder.
#
# This file does only contain a selection of the most common options. For a
# full list see the documentation:
# http://www.sphinx-doc.org/en/master/config

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
# import os
# import sys
# sys.path.insert(0, os.path.abspath('.'))

# -- Import dependencies -----------------------------------------------------

import sphinx_rtd_theme # the readthedocs theme module

# -- Customize code highlighting ---------------------------------------------

from pygments.style import Style
from pygments.styles import get_all_styles
from pygments.styles.pastie import PastieStyle
from pygments.lexer import inherit, bygroups, include, combined
from pygments.lexers.c_cpp import CppLexer
from pygments.lexers.python import PythonLexer, Python3Lexer
from pygments.token import *
from sphinx.highlighting import lexers

class ExtendedCppLexer(CppLexer):
    tokens = {
        'statements': [
            include('xfunction'),
            include('xmethod'),
            include('xclass'),
            inherit
        ],
        'xfunction': [
            (r'(\s*)([^\w\s\d\"\']+)(\s*)([a-zA-Z][a-zA-Z0-9_]*)(\s*)(\()',
                bygroups(Whitespace, Operator, Whitespace, Name.Function, Whitespace, Punctuation)),
        ],
        'xmethod': [
            (r'([a-zA-Z][a-zA-Z0-9_]*)(\s*)(\.)(\s*)([a-zA-Z][a-zA-Z0-9_]*)(\s*)(\()',
                bygroups(Name, Whitespace, Punctuation, Whitespace, Name.Function, Whitespace, Punctuation)),
        ],
        'xclass': [
            (r'([A-Z][a-zA-Z0-9_]*)(\s+)([a-zA-Z0-9_]+)',
                bygroups(Name.Class, Whitespace, Name)),
        ],
    }

class ExtendedPythonLexer(PythonLexer):
    tokens = {
        'name': [
            include('xfunction'),
            include('xmethod'),
            include('xclass'),
            inherit
        ],
        'xfunction': [
            (r'(\s*)([a-z][a-zA-Z0-9_]*)(\s*)(\()',
                bygroups(Whitespace, Name.Function, Whitespace, Punctuation)),
        ],
        'xmethod': [
            (r'([a-zA-Z][a-zA-Z0-9_]*)(\s*)(\.)(\s*)([a-zA-Z][a-zA-Z0-9_]*)(\s*)(\()',
                bygroups(Name, Whitespace, Punctuation, Whitespace, Name.Function, Whitespace, Punctuation)),
        ],
        'xclass': [
            (r'([a-zA-Z0-9_]+)(\s*)(\=)(\s*)([A-Z][a-zA-Z0-9_]+)',
                bygroups(Name, Whitespace, Operator, Whitespace, Name.Class)),
        ],
    }

class ExtendedPython3Lexer(ExtendedPythonLexer, Python3Lexer):
   pass

# Set the extended and customized C++ and Python lexers
lexers['cpp'] = lexers['c++'] = ExtendedCppLexer()
lexers['python'] = lexers['py'] = ExtendedPythonLexer()
lexers['python3'] = lexers['py3'] = ExtendedPython3Lexer()

class ReaktoroStyle(PastieStyle):
    background_color = '#F8F8F8'

    highlight_color = '#FFD7BC'

    default_style = ''

    styles = {

        Whitespace:             '#bbbbbb',
        Comment:                '#888888',
        Comment.Preproc:        'bold #e74c3c', # 'bold #cc0000',
        Comment.Special:        'bg:#fff0f0 bold #e74c3c', # 'bg:#fff0f0 bold #cc0000',

        String:                 '#9b59b6', # 'bg:#fff0f0 #dd2200',
        String.Regex:           '#2c3e50', # 'bg:#fff0ff #008800',
        String.Other:           '#22bb22', # 'bg:#f0fff0 #22bb22',
        String.Symbol:          '#aa6600',
        String.Interpol:        '#3333bb',
        String.Escape:          '#0044dd',

        Operator.Word:          '#2c3e50',

        Keyword:                'bold #2c3e50',
        Keyword.Pseudo:         'nobold',
        Keyword.Type:           '#888888',

        Name.Class:             'bold #e74c3c', # 'bold #bb0066',
        Name.Exception:         'bold #e74c3c', # 'bold #bb0066',
        Name.Function:          'bold #2980b9', # 'bold #0066bb',
        Name.Property:          'bold #336699',
        Name.Namespace:         'bold #e74c3c', # 'bold #bb0066',
        Name.Builtin:           '#003388',
        Name.Variable:          '#336699',
        Name.Variable.Class:    '#336699',
        Name.Variable.Instance: '#3333bb',
        Name.Variable.Global:   '#dd7700',
        Name.Constant:          'bold #003366',
        Name.Tag:               'bold #e74c3c', # 'bold #bb0066',
        Name.Attribute:         '#336699',
        Name.Decorator:         '#555555',
        Name.Label:             'italic #336699',

        Number:                 '', # 'bold #0000DD',

        Generic.Heading:        '#333',
        Generic.Subheading:     '#666',
        Generic.Deleted:        'bg:#ffdddd #000000',
        Generic.Inserted:       'bg:#ddffdd #000000',
        Generic.Error:          '#aa0000',
        Generic.Emph:           'italic',
        Generic.Strong:         'bold',
        Generic.Prompt:         '#555555',
        Generic.Output:         '#888888',
        Generic.Traceback:      '#aa0000',

        Error:                  'bg:#e3d2d2 #a61717'
    }

# Create on-the-fly a reaktoro module under
# pygments.styles with the ReaktoroStyle class.
# Reference: https://stackoverflow.com/questions/48615629/how-to-include-pygments-styles-in-a-sphinx-project?rq=1
import sys
import pygments.styles
sys.modules['pygments.styles.reaktoro'] = type(sys)('reaktoro')
sys.modules['pygments.styles.reaktoro'].ReaktoroStyle = ReaktoroStyle
pygments.styles.reaktoro = 'ReaktoroStyle'
pygments.styles.STYLE_MAP['reaktoro'] = 'reaktoro::ReaktoroStyle'

# -- Project information -----------------------------------------------------

project = 'Reaktoro'
copyright = '2021, Allan Leal and Reaktoro contributors'
author = 'Allan Leal'

# The short X.Y version
version = ''
# The full version, including alpha/beta/rc tags
release = ''


# -- General configuration ---------------------------------------------------

# If your documentation needs a minimal Sphinx version, state it here.
#
# needs_sphinx = '1.0'

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'sphinx.ext.autosectionlabel',
    'sphinx.ext.mathjax',
    'sphinx.ext.todo',
    'sphinxcontrib.images',
]

# If this is True, todo and todolist produce output, else they produce nothing.
todo_include_todos = True

# If this is True, todo emits a warning for each TODO entries.
todo_emit_warnings = True

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# The suffix(es) of source filenames.
# You can specify multiple suffix as a list of string:
#
# source_suffix = ['.rst', '.md']
source_suffix = '.rst'

# The master toctree document.
master_doc = 'index'

# The language for content autogenerated by Sphinx. Refer to documentation
# for a list of supported languages.
#
# This is also used if you do content translation via gettext catalogs.
# Usually you set "language" from the command line for these cases.
language = None

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

# The name of the Pygments (syntax highlighting) style to use.
pygments_style = 'reaktoro'

# -- reST epilog and prolog Options for HTML output --------------------------
rst_epilog = """
"""
rst_prolog = """
"""

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
html_theme = 'sphinx_rtd_theme'
html_theme_path = [sphinx_rtd_theme.get_html_theme_path()]

# Theme options are theme-specific and customize the look and feel of a theme
# further.  For a list of options available for each theme, see the
# documentation.

html_theme_options = {
    'titles_only': True,
}

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

# Custom sidebar templates, must be a dictionary that maps document names
# to template names.
#
# The default sidebars (for documents that don't match any pattern) are
# defined by theme itself.  Builtin themes are using these templates by
# default: ``['localtoc.html', 'relations.html', 'sourcelink.html',
# 'searchbox.html']``.
#
# html_sidebars = {}


# -- Options for HTMLHelp output ---------------------------------------------

# Output file base name for HTML help builder.
htmlhelp_basename = 'Reaktorodoc'


# -- Options for LaTeX output ------------------------------------------------

latex_elements = {
    # The paper size ('letterpaper' or 'a4paper').
    #
    # 'papersize': 'letterpaper',

    # The font size ('10pt', '11pt' or '12pt').
    #
    # 'pointsize': '10pt',

    # Additional stuff for the LaTeX preamble.
    #
    # 'preamble': '',

    # Latex figure (float) alignment
    #
    # 'figure_align': 'htbp',
}

# Grouping the document tree into LaTeX files. List of tuples
# (source start file, target name, title,
#  author, documentclass [howto, manual, or own class]).
latex_documents = [
    (master_doc, 'Reaktoro.tex', 'Reaktoro Documentation',
     'Allan Leal', 'manual'),
]


# -- Options for manual page output ------------------------------------------

# One entry per manual page. List of tuples
# (source start file, name, description, authors, manual section).
man_pages = [
    (master_doc, 'reaktoro', 'Reaktoro Documentation',
     [author], 1)
]


# -- Options for Texinfo output ----------------------------------------------

# Grouping the document tree into Texinfo files. List of tuples
# (source start file, target name, title, author,
#  dir menu entry, description, category)
texinfo_documents = [
    (master_doc, 'Reaktoro', 'Reaktoro Documentation',
     author, 'Reaktoro', 'One line description of project.',
     'Miscellaneous'),
]


# -- Options for Epub output -------------------------------------------------

# Bibliographic Dublin Core info.
epub_title = project

# The unique identifier of the text. This can be a ISBN number
# or the project homepage.
#
# epub_identifier = ''

# A unique identification for the text.
#
# epub_uid = ''

# A list of files that should not be packed into the epub file.
epub_exclude_files = ['search.html']


# -- Extension configuration -------------------------------------------------
