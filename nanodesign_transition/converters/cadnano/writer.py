#!/usr/bin/env python
"""
This module is used to write caDNAno design JSON files.
"""
import json
import logging
from cadnano.common import CadnanoLatticeName,CadnanoLatticeType

try:
    import os.path
    import sys
    base_path = os.path.abspath( os.path.dirname(__file__) + '/../../../' )
    sys.path.append( base_path )
    from nanodesign_transition.parameters import DnaParameters
    sys.path = sys.path[:-1]
except ImportError as i:
    print "Could not get nanodesign_transition module"
    raise i

class CadnanoWriter(object):
    """ The CadnanoWriter class is used to write out a caDNAno design JSON file.
    """
    def __init__(self, dna_structure):
        self.dna_structure = dna_structure
        self._logging_level = logging.INFO
        self._setup_logging()

    def _setup_logging(self):
        """ Set up logging."""
        self._logger = logging.getLogger(__name__)
        #self._logger = logging.getLogger('viewer:writer')
        self._logger.setLevel(self._logging_level)

        # Create console handler and set format.
        console_handler = logging.StreamHandler()
        #formatter = logging.Formatter('%(asctime)s [%(name)s] %(levelname)s - %(message)s')
        formatter = logging.Formatter('[%(name)s] %(levelname)s - %(message)s')
        console_handler.setFormatter(formatter)
        self._logger.addHandler(console_handler)

    def write(self,file_name):
        """ Write a caDNAno design JSON file.

        Args:
            file_name (string): The name of a caDNAno JSON file to write.
        """
        self._logger.info("Writing caDNAno design JSON file: %s " % file_name)
        dna_structure = self.dna_structure
        dna_structure.compute_aux_data()
        vstrand_info = self._get_vstrand_info(dna_structure)

        design = { 'name'    : os.path.basename(file_name),
                   'vstrands': vstrand_info
                 }

        with open(file_name, 'w') as outfile:
            json.dump(design, outfile, indent=4, separators=(',', ': '))

    def _get_vstrand_info(self, dna_structure):
        """ Get virtual helix information for the design. 
            In caDNAno all the virtual helices are the same size (i.e. the same number of base positions).
            The lists of bases for staple and scaffold strands are also the same size. Positions within
            this list that contain no bases are set to None. 
        """
        self._logger.info("Number of helices %d " % len(dna_structure.structure_helices_map))
        helix_size = 0
        vstrands_info = []

        # Create a map of the helices by ID so we can reproduce the order in which 
        # they were read in from the original caDNAno file.
        helix_map = {}
        for helix in dna_structure.structure_helices_map.itervalues():
            helix_map[helix.id] = helix 

        # Get the base connectivity, loops and skips for the staple and scaffold strands.
        for id in sorted(helix_map):
            helix = helix_map[id] 
            scaffold_base_list = helix.scaffold_base_list
            staple_base_list = helix.staple_base_list
            helix_size = len(scaffold_base_list)
            id = helix.id 
            row = helix.lattice_row 
            col = helix.lattice_col
            num = helix.lattice_num
            self._logger.info("Helix id %d  num %d  row %d  col %d  size %d" % (id, num, row, col, helix_size)) 

            scaf_info = self._get_base_info(scaffold_base_list)
            stap_info = self._get_base_info(staple_base_list)
            loop_info = self._get_loop_info(staple_base_list)
            skip_info = self._get_skip_info(staple_base_list)

            vstrand = { "row": row,
                        "col": col,
                        "num": num,
                        "scaf": scaf_info,
                        "stap": stap_info,
                        "loop": loop_info,
                        "skip": skip_info,
                        "stap_colors": []
                      }

            vstrands_info.append( vstrand )

        return vstrands_info

    def _get_base_info(self, base_list):
        """ Get the base connecivity information from the list of bases. """
        base_info = []
        for base in base_list:
            if base:
                #print(">>> base  id %d   h %d  p %d  up %d  down %d " % (base.id, base.h, base.p, base.up, base.down))
                if base.up == None:
                    up_pos = -1
                    up_vh = -1
                else:
                    up_pos = base.up.p
                    up_vh = base.up.h

                if base.down == None:
                    down_pos = -1
                    down_vh = -1
                else:
                    down_pos = base.down.p
                    down_vh = base.down.h
                #print("          %s " % str([up_vh, up_pos, down_vh, down_pos])) 
                base_info.append( [up_vh, up_pos, down_vh, down_pos] )
            else:
                base_info.append( [-1,-1,-1,-1] )
        #__for base in base_list
        return base_info 

    def _get_loop_info(self, base_list):
        """ Get the loop (insert) information for each base. """
        loop_info = []
        for base in base_list:
            if base:
                loop_info.append(base.num_insertions)
            else:
                loop_info.append(0)
        return loop_info

    def _get_skip_info(self, base_list):
        """ Get the skip (deletion) information for each base. """
        skip_info = []
        for base in base_list:
            if base:
                skip_info.append(base.num_deletions)
            else:
                skip_info.append(0)
        return skip_info
