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
caDNAno common definitions.

This module defines some common definitions used to process caDNAno DNA origami design JSON files. 
"""

class CadnanoStrandType:
    SCAFFOLD = 0
    STAPLE   = 1

class CadnanoLatticeName:
    """Cadnano lattice names"""
    NONE      = "none"
    HONEYCOMB = "honeycomb"
    SQUARE    = "square"

class CadnanoLatticeType:
    """Cadnano lattice types"""
    none      = 0
    honeycomb = 1
    square    = 2
    names = [CadnanoLatticeName.NONE, CadnanoLatticeName.HONEYCOMB, CadnanoLatticeName.SQUARE]

class CadnanoJsonFields:
    """Field names in caDNAno json file.""" 
    COL         = "col"
    DELETION    = "deletion"
    INSERTION   = "insertion"
    LOOP        = "loop"
    NUM         = "num"
    ROW         = "row"
    SCAF        = "scaf"
    SCAFFOLD    = "scaffold"
    SKIP        = "skip"
    STAPLE      = "staple"
    STAP        = "stap"
    STAP_COLORS = "stap_colors"
    VHELIX      = "vhelix"
    VSTRANDS    = "vstrands"

