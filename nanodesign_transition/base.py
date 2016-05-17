"""
This module is used to store data for a DNA base. 

A list of DnaBase objects, stored in the DnaStructure object, is used to define the connectivity (topology) 
of a DNA structure: bases that are adjacent along a single DNA helix and which base they are paired with (if any). 
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
            across (int): The base ID of the base's Watson-Crick neighbor.
            domain (int): The domain ID the base is in.
            down (int): The base ID of the base's 3' neighbor.
            h (int): The ID of the helix the base is in. 
            id (int): The base ID.
            is_scaf (bool): If True then this base is in a scaffold strand. 
            loop (int): The number of insertions at this base.
            p (int): The helix position of the base.
            seq (string): A one character string representing the base sequence nucleotide. 
            skip (int): The number of deletions at this base.
            strand (int): The strand ID the base is in.
            up (int): The base ID of the base's 5' neighbor.
    """

    def __init__( self, id, up=0, down=0, across=0, seq='N'):
        self.id = int(id)
        self.up = int(up)
        self.down = int(down)
        self.across = int(across)
        self.seq = seq
        self.strand = -1
        self.domain = -1
        self.loop = 0
        self.skip = 0
        self.h = -1
        self.p = -1
        self.is_scaf = False

