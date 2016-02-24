"""
NanoDesign
==========

Provides
  1. A interface into multiple file types found in various DNA nanostructure design tools, such as Cadnano, (others?)
  2. A data format for representing DNA nanostructure designs, suitable for visualization and algorithmic use.
  3. Standard algorithms to run on DNA nanostructure designs, for characterizing various statistics or performing editing operations.

Documentation is currently sparse at best, but you should be able to explore the docstrings of various functions and packages.

Available subpackages
---------------------

data
    Internal data formats
converters
    File loading and conversion routines.
algorithms
    Statistics and other algorithms on these file types (maybe should just be in the core package? maybe all of it should?)

"""
from __future__ import print_function

# Load the basic core elements, as well as subpackages.
from . import core
from .core import *

from . import data
from . import converters
from . import algorithms

# Designate which components will be in the * namespace.
__all__ = []
__all__.extend(core.__all__)
__all__.extend(['data','converters','algorithms'])

