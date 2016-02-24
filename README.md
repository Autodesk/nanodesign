NanoDesign
==========

This repository contains the NanoDesign python package, as well as supporting material, such as docs, tests, Dockerfile, scripts used for common tasks, scripts used for the NanoDesign frontend tasks, etc.

## Authors
 
For more information, contact:

* Joseph Schaeffer <joseph.schaeffer@autodesk.com>
* Dave Parker <dave.parker@autodesk.com>
* Mike Zyracki <michael.zyracki@autodesk.com>

## Directory Structure

We use the following directory structure:

* `nanodesign`: the main python package
* `docs`: doc generation for the package
* `tests`: tests for the package
* `scripts`: example scripts for using the package, as well as scripts used to run specific tasks for the NanoDesign web application

Individual top level files of interest:

* `README.md`: You're reading it.
* `Makefile`: Some basic commands you might use, such as 'make init' for doing the requirements install, or 'make tests' for running the tests. More will be added.
* `requirements.txt`: pip requirements file. Add specific package requirements here.
* `setup.py`: distutils based package install. To run the package install, 'python setup.py install'
* `LICENSE`: Our distribution license. Currently blank, do not distribute.
* `Dockerfile`: Defines docker container with appropriate requirements.

## Purpose

### Python Package
 We want to build most of our utilities as a python package; this will enable later API releases and open-source coding to expand the types of conversions allowed. The package is envisioned as having several basic components:

* Data Format(s) - our internal data format for use by specific algorithms and to serve as an intermediate if needing to do format conversion. Could be multiples depending on usage. Should also have the format for sending back to client, etc.
* Statistics Algorithms - python implementations of various visualization components being discussed on the frontend. This would let us make those easily expandable (if slower) by allowing others to be added on the backend if our data format can support sending it.
* Editing/Modification Algorithms - e.g. auto-stapling, perhaps specific shape auto-creation, etc.
* Converters - allow file type conversion and input to our system. Cadnano, CanDo, etc.


