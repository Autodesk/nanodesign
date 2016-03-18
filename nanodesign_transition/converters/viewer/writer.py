#!/usr/bin/env python
""" 
This module is used to write DNA Design viewer JSON files. 
"""
import json
import logging
import numpy as np
from cadnano.common import CadnanoLatticeName,CadnanoLatticeType
from nanodesign.parameters import DnaParameters

class ViewerWriter(object):
    """ The ViewerWriter class writes out a DNA Design viewer JSON file. 
    """

    def __init__(self, dna_structure):
        self.dna_structure = dna_structure
        self._logging_level = logging.INFO
        self._setup_logging()

    def set_logging_level(self,level):
        """Set logging level."""
        self._logger.setLevel(level)

    def write(self,file_name):
        """Write a viewer JSON file.

        Args:
            file_name (string): The name of a viewer JSON file to write.

        """

        self._logger.info("Writing DNA Design Viewer JSON file: %s " % file_name)
        dna_structure = self.dna_structure
        dna_structure.compute_aux_data()
        lattice_type = CadnanoLatticeType.names[dna_structure.lattice_type]
        strands_info = self._get_strand_info(dna_structure)
        helices_info = self._get_helices_info(dna_structure)
        domains_info = self._get_domain_info(dna_structure)

        vis_model = { 'model_name'      : dna_structure.name,
                      'lattice_type'    : lattice_type,
                      'strands'         : strands_info,
                      'virtual_helices' : helices_info,
                      'domains'         : domains_info
                    }

        with open(file_name, 'w') as outfile:
            json.dump(vis_model, outfile, indent=4, separators=(',', ': '))
            #json.dump(vis_model, outfile)

    def _get_domain_info(self, dna_structure):
        """ Get JSON serialized data for all the domains. """
        domains_info = [] 
        for domain in dna_structure.domain_list:
            point1,point2 = domain.get_end_points()
            base_info = [base.id for base in domain.base_list]

            if (domain.strand):
                strand_id = domain.strand.id
                strand_bases = domain.strand.tour
                start_base_index = strand_bases.index(base_info[0])
                end_base_index = strand_bases.index(base_info[-1])
            else:
                strand_id = -1
                start_base_index = -1
                end_base_index = -1

            frame = domain.helix.end_frames[:,:,0]

            info = { 'id'                : domain.id ,
                     'color'             : domain.color,
                     'strand_id'         : domain.strand.id,
                     'start_position'    : [point1[0],point1[1],point1[2]],
                     'end_position'      : [point2[0],point2[1],point2[2]],
                     'orientation'       : [frame[0,2], frame[1,2], frame[2,2]],
                     'number_of_bases'   : len(domain.base_list),
                     'bases'             : base_info,
                     'start_base_index'  : start_base_index,
                     'end_base_index'    : end_base_index,
                     'connected_strand'  : domain.connected_strand,
                     'connected_domain'  : domain.connected_domain
                   }
            domains_info.append(info)
        #__for domain in dna_structure.domain_list
        return domains_info 

    def _get_helices_info(self, dna_structure):
        """ Get JSON serialized data for helix objects. """
        helices_info = []
   
        for helix in dna_structure.structure_helices:
            point1 = helix.end_coordinates[0]
            point2 = helix.end_coordinates[1]
            length = np.linalg.norm(point1-point2)
            frame = helix.end_frames[:,:,0]
            domain_info = [ domain.id for domain in helix.domain_list ]

            info = { 'id'                : helix.id,
                     'length'            : length,
                     'strand_radius'     : DnaParameters.strand_radius,
                     'base_pair_rise'    : DnaParameters.base_pair_rise,
                     'start_position'    : list(point1),
                     'orientation'       : [frame[0,2], frame[1,2], frame[2,2]],
                     'end_position'      : list(point2),
                     'scaffold_polarity' : helix.scaffold_polarity,
                     'cadnano_info'      : { 'row' : helix.lattice_row, 'col' : helix.lattice_col, 
                                             'num' : helix.lattice_num },
                     'domains' : domain_info
                   }
            helices_info.append(info)
        return helices_info

    def _get_strand_info(self, dna_structure):
        """ Get JSON serialized data for strands objetcs. 
        """
        strand_info_list = []
        for strand in dna_structure.strands:
            domains_info = strand.get_domains_info()
            base_coords = strand.get_base_coords()
            base_info = []
            for i in xrange(0,len(strand.tour)):
                id = strand.tour[i]
                base = self.dna_structure.base_connectivity[id-1]
                coord = base_coords[i]
                base_info.append({ 'id': base.id, 'coordinates' : list(coord), 'sequence' : base.seq })
 
            info = { 'id' : strand.id,
                     'is_scaffold'     : strand.is_scaffold,
                     'is_circular'     : strand.is_circular,
                     'number_of_bases' : len(strand.tour),
                     'virtual_helices' : list(strand.helix_list.keys()),
                     'domains'         : domains_info,
                     'bases'           : base_info,
                     'color'           : strand.color
                   }

            strand_info_list.append(info)
        return strand_info_list

    def _setup_logging(self):
        """ Set up logging."""
        self._logger = logging.getLogger('viewer:writer')
        self._logger.setLevel(self._logging_level)

        # create console handler and set format
        console_handler = logging.StreamHandler()
        #formatter = logging.Formatter('%(asctime)s [%(name)s] %(levelname)s - %(message)s')
        formatter = logging.Formatter('[%(name)s] %(levelname)s - %(message)s')
        console_handler.setFormatter(formatter)
        self._logger.addHandler(console_handler)

def main():
    """ Write a DNA Design Viewer JSON file."""
    file_name = sys.argv[1]
    cadnano_reader = CadnanoReader()
    #cadnano_reader.set_logging_level(logging.DEBUG)
    cadnano_reader.read_json(file_name)

if __name__ == '__main__':
    main()

