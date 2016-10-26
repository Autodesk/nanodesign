#!/usr/bin/env python
"""
This module is used to store information for a DNA strand.

A DNA strand is a continuous chain of nucleotides. It can be either a scaffold of a staple.
"""
import sys
import os
from sets import Set
import numpy as np

class DnaStrand(object):
    """ The DnaStrand class stores data for a DNA strand.

    Attributes:
        base_id_list (Dict): The location of each strand base within list of base IDs making up the strand. 
            The dictionary maps base IDs to an index into tour[].
        color (List[float]: The strand color in RGB.
        dna_structure (DnaStructure): The DNA structure this strand belongs to. 
        domain_list (List[Domain]): The list of domains for this strand.
        helix_list (Dict): The list of helices the strand passes through. The dictionary maps helix IDs 
            to DnaStructureHelix objects.
        icolor (int): The strand color as an integer. The integer color can be used as an ID to group staple strands. 
        id (int): The strand ID. 
        insert_seq (List[string]): The list of sequence letters inserted into this strand. 
        is_circular (bool): If True then the strand is circular, returning to its starting postion.
        is_scaffold (bool): If True then the strand is a scaffold strand.
        tour (List[int]): The list of base IDs making up the strand. 
    """

    def __init__(self, id, dna_structure):
        self.id = id
        self.is_scaffold = False
        self.is_circular = False
        self.tour = []
        self.color = [1.0,1.0,1.0]
        self.icolor = None
        self.helix_list = dict()
        self.base_id_list = dict()
        self.dna_structure = dna_structure
        self.domain_list = []
        self.insert_seq = []

    def add_helix(self, helix): 
        id = helix.lattice_num
        if (id not in self.helix_list):
            #print("[DnaStrand] ---------- strand %d  add helix %d ----------" % (self.id, id))
            self.helix_list[id] = helix
    #__def add_helix

    def get_base_coords(self):
        """ Get the coordinates of bases along the dna helix axis. 
            This is only used when writing a visualization file.
        """
        num_bases = len(self.tour)
        base_coords = np.zeros((num_bases,3), dtype=float)
        for i,base in enumerate(self.tour):
            helix_num = base.h
            helix_pos = base.p
            helix = self.helix_list[helix_num]
            base_coords[i] = base.coordinates
        return base_coords
    #__def get_base_coords

    def get_base_index(self, base):
        """ Get the index into the strand for the given base. 
            This is only used when writing a visualization file.
        """
        num_bases = len(self.tour)
        if (not self.base_id_list):
            for i,sbase in enumerate(self.tour):
                self.base_id_list[sbase.id] = i
        if base.id not in self.base_id_list:
            sys.stderr.write("[strand::get_base_index] **** WARNING: base %d not found in strand %d.\n" % (base.id, self.id))
            return None
        return self.base_id_list[base.id]
    #__def get_base_index

