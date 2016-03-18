#!/usr/bin/env python
"""
This module is used to store data for a DNA base. 

A list of DnaBase objects, stored in the DnaStructure object, is used to define the connectivity (topology) 
of a DNA structure: bases that are adjacent along a single DNA helix and whixh base they are paired with (if any). 
"""
import os
import sys
import logging
import re
import string
import sys
import time
import numpy as np

class DnaBase(object):
    """ This class stores information for a DNA base. 

        Attributes:
            id (int): ID
            up (int): 5' neighbor
            down (int): 3' neighbor
            across (int): The Watson-Crick neighbor base ID.
            coord (numpy.ndarray shape=(3,)): 3D coordinate on the DNA helix backbone 
            skip (int): deletion
            loop (int): insertion
            h (int): Helix ID (start from 0)
            p (int): Position in a helix (start from 0)
            is_scaf (bool): If this base is on a scaffold. 
            strand (int): 
            residue (int):
            seq (string): A one character string representing the base sequence nucleotide. 
    """

    def __init__( self, id, up=0, down=0, across=0, seq='N'):
        self.id = int(id)
        self.up = int(up)
        self.down = int(down)
        self.across = int(across)
        self.seq = seq
        self.strand = -1
        self.domain = -1
        self.residue = 0
        self.loop = 0
        self.skip = 0
        self.h = -1
        self.p = -1
        self.is_scaf = False
        self.coord = np.array([0.0,0.0,0.0],dtype=float)

    def set_domain_id(self, id):
        self.domain = id

