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
This module is used to store data read from caDNAno DNA origami design JSON files. 
"""
import logging
import sys
from itertools import product

from .common import CadnanoLatticeType,CadnanoStrandType
from ...data.base import DnaBase
from ...data.strand import DnaStrand
from ...data.dna_structure import DnaStructure,DnaStructureHelix
from ...data.lattice import Lattice,SquareLattice,HoneycombLattice


class CadnanoDesign(object):
    """The CadnanoDesign class stores data read from a caDNAno origami design JSON file."""

    def __init__(self, name="design"):
        self.name = name
        self.max_base_id = 0
        self.lattice_type = CadnanoLatticeType.none
        self.helices = []
        self.helices_coord_map = {}
        self.max_row = 0
        self.max_col = 0
        self._logger = logging.getLogger(__name__) 

    def calculate_possible_crossovers(self):
        """ Calculate the possible crossover positions for both scaffold and staples. 

            The method for computing cross-ovrs is based on that from the caDNAno part.py file. 
        """
        self._logger.debug(" ======================== calculate crossovers ======================== ")
        self._logger.debug(">>> lattice type: %s " % CadnanoLatticeType.names[self.lattice_type])

        if self.lattice_type == CadnanoLatticeType.square:
            lattice = SquareLattice 
        else:
            lattice = HoneycombLattice 

        strand_types = (CadnanoStrandType.SCAFFOLD, CadnanoStrandType.STAPLE)

        luts_neighbor = list( zip( lattice.scaffold_low,
                                   lattice.scaffold_high,
                                   lattice.staple_low ,
                                   lattice.staple_high ))

        num_bases = self.max_base_id
        base_range = list(range(0, num_bases, lattice.step))
        self._logger.debug(">>> num_bases: %d" % num_bases)
        self._logger.debug(">>> base_range_unit: %s" % str(base_range))
        self._logger.debug(">>> luts_neighbor: %s" % str(luts_neighbor))

        # Create a mapping from lattice coordinates to helix.
        for vhelix in self.helices:
            self.helices_coord_map[(vhelix.row,vhelix.col)] = vhelix

        for vhelix in self.helices:
            self._logger.debug("")
            self._logger.debug("-------------------------- helix num %d -------------------------" % vhelix.num)

            neighbor_helices = self.get_neighbor_helices(lattice, vhelix)
            for neighbor, lut in zip(neighbor_helices, luts_neighbor):
                if not neighbor:
                    continue
                self._logger.debug("----------- neighbor %s -----------" % str(neighbor.num))
                lut_scaf = lut[0:2]
                lut_stap = lut[2:4]
                lut = (lut_scaf, lut_stap)
                self._logger.debug(">>> lut_stap: %s" % str(lut_stap))

                for pts, st in zip(lut, strand_types):
                    if st == CadnanoStrandType.SCAFFOLD:
                        stype = "scafflold"
                    else:
                        stype = "staple"

                    for pt in pts:
                        self._logger.debug(">>>>>> pt: %s   st: %s" % (str(pt),stype))
                        self._logger.debug(">>>>>> base_range: %s " % (str(base_range)))
                        for i, j in product(base_range, pt):
                            index = i + j
                            if index < num_bases:
                                if st == CadnanoStrandType.SCAFFOLD:
                                    vhelix.possible_scaffold_crossovers.append((neighbor,index))
                                else:
                                    vhelix.possible_staple_crossovers.append((neighbor,index))
                                self._logger.debug("    add i: %4d  j: %4d  index: %4d" % (i,j,index))
                    #__for pt in pts
                #__for pts, st in zip(lut, strand_types)
            #__for neighbor, lut in zip(neighbor_helices, luts_neighbor):
        #__for vhelix in self.virtual_helices

        #sys.exit(0)

    def get_neighbor_helices(self, lattice, vhelix):
        """ Get the neighboring helices for the given helix.
        """
        neighbor_helices = []
        row = vhelix.row
        col = vhelix.col
        neighbor_coords = lattice.get_neighbors(row,col)
        for coord in neighbor_coords:
            if coord in self.helices_coord_map:
                nhelix = self.helices_coord_map[coord]
                neighbor_helices.append(nhelix)
            else:
                neighbor_helices.append(None)
        #__for coord in neighbor_coords
        return neighbor_helices


class CadnanoVirtualHelix(object):
    """ This class stores data for a caDNAno virtual helix. 

        Attributes: 
            id (int): The id of the virtual helix.
            num (int): The index of the virtual helix. This is the row number in the caDNAno 2D design diagram.
            col (int): The column index of the virtual helix in a lattice. 
            row (int): The row index of the virtual helix in a lattice. 
            possible_staple_crossovers (list[(CadnanoVirtualHelix,int)]: The list of possible staple crossovers to 
                neighboring helices. 
            possible_scaffold_crossovers (list[(CadnanoVirtualHelix,int)]: The list of possible scaffold crossovers to 
                neighboring helices.  
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
        self.possible_staple_crossovers = []
        self.possible_scaffold_crossovers = []


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

