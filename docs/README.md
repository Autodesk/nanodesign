Documentation
=============

Most of the entries in this directory were generated via sphinx-quickstart. Relevant files for editing:

* `conf.py`: the generated configuration for the sphinx documentation, may be modified if we start changing settings, etc.
* `index.rst`: the main source file for the documentation.

Probably shouldn't edit:
* `Makefile`: has the basic commands for building the docs. E.g. for an html version, 'make html' will create it.
* `make.bat`: same as Makefile, for Windows users.


Usage 
-----

Make sure sphinx is installed. `make init` in the root directory should install it, but if not, just use `pip install sphinx`. 
Once sphinx is installed, `make html` in this directory will build the html version of the documentation. It will then reside in `_build/html/index.html` as the root. 