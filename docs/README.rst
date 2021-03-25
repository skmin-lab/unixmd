***************************************
Compile the PyUNIxMD Manual with Sphinx
***************************************

Requirements
============
* Sphinx
* sphinx_rtd_theme
* sphinxcontrib-bibtex

if you don't have Sphinx,

::

  $ pip install -U Sphinx

if you don't have sphinx_rtd_theme,

::

  $ pip install sphinx_rtd_theme

if you don't have sphinxcontrib-bibtex,

::

  $ pip install sphinxcontrib-bibtex

Build
=====    
You can locally build the manual by the following command.

::

  $ make html

Test
====    
You can check out the manual by the following command.

::

  $ cd build/html
  $ firefox unixmd.html

