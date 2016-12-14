# Copyright 2016 Autodesk Inc.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

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
from . import utils
from . import visualizer

from .data import Domain

from .data import energy_model
from .data import convert_temperature_K_to_C

# Create a logger console handler and set logging output format.
def _init_logging():
    import logging
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.INFO)
    console_handler = logging.StreamHandler()
    formatter = logging.Formatter('[%(name)s] %(levelname)s - %(message)s')
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)

_init_logging()

# Designate which components will be in the * namespace.
__all__ = []
__all__.extend(core.__all__)
__all__.extend(['data','converters','algorithms','utils','visualizer'])

