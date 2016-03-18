#!/usr/bin/env python
""" 
caDNAno common definitions.

This module defines some common definitions used to process caDNAno DNA origami design JSON files. 
"""

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

