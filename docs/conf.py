import os
project = 'GRABNES'
author = 'Jeil Jung Group'
extensions = [
    'breathe',
    'myst_parser',
    'sphinx.ext.mathjax',
]

templates_path = ['_templates']
exclude_patterns = [
    '_build',
    'development/pre-announcement-changelog.md',
]
html_theme = 'sphinx_rtd_theme'

# Breathe configuration
breathe_projects = {
    'grabnes': os.path.abspath('_build/doxygen/xml'),
}
breathe_default_project = 'grabnes'
