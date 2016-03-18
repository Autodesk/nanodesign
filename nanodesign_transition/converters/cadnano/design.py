#!/usr/bin/env python
""" 
This module is used to store data read from caDNAno DNA origami design JSON files. 
"""
import logging
from common import CadnanoLatticeType

class CadnanoDesign(object):
    """The CadnanoDesign class stores data read from a caDNAno origami design JSON file."""

    def __init__(self, name="design"):
        self.name = name
        self.num_bases = 0
        self.lattice_type = CadnanoLatticeType.none
        self.helices = []
        self.max_row = 0
        self.max_col = 0
        self._logging_level = logging.INFO
        self._logging = self._setup_logging()

    def set_logging_level(self,level):
        """Set logging level."""
        self._logger.setLevel(level)

    def _setup_logging(self):
        """ Set up logging."""
        self._logger = logging.getLogger('cadnano:model')
        self._logger.setLevel(self._logging_level)

        # create console handler and set format
        console_handler = logging.StreamHandler()
        #formatter = logging.Formatter('%(asctime)s [%(name)s] %(levelname)s - %(message)s')
        formatter = logging.Formatter('[%(name)s] %(levelname)s - %(message)s')
        console_handler.setFormatter(formatter)
        self._logger.addHandler(console_handler)

class CadnanoVirtualHelix(object):
    """ This class stores data for a caDNAno virtual helix. 

        Attributes: 
            id (int): The id of the virtual helix.
            num (int): The index of the virtual helix. This is the row number in the caDNAno 2D design diagram.
            col (int): The column index of the virtual helix in a lattice. 
            row (int): The row index of the virtual helix in a lattice. 
    """
    def __init__(self, id, num, row, col, insertions, deletions):
        self.id = id
        self.num = num
        self.row = row
        self.col = col
        self.insertions = insertions
        self.deletions = deletions
        self.scaffold_strands = []
        self.staple_strands = []
        self.staple_colors = []

class CadnanoBase():
    """ This class stores data for a caDNAno base, either staple or scaffold. 

        The data here is equivalent to each 4-tuple [V_0,b_0,V_1,b_1] entry stored in the caDNAno JSON
        'scaf' and 'stap' arrays.

        Attributes: 
            initial_strand (int): id of the initial virtual helix the base connects to. 
            initial_base (int): id of the initial base.
            final_strand (int): id of the final virtual helix the base connects to. 
            final_base (int): id of the final base. 
        
    """
    def __init__(self, initial_strand, initial_base, final_strand, final_base):
        self.initial_strand = initial_strand
        self.initial_base = initial_base
        self.final_strand = final_strand
        self.final_base = final_base 

