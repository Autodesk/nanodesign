#!/usr/bin/env python
""" 
This module is used to write DNA Design viewer JSON files. 
"""
import collections
import json
import logging
import numpy as np
from math import pi
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
                     'connected_domain'  : domain.connected_domain,
                     'melting_temperature' : domain.melting_temperature() #"{:.2f}".format(domain.melting_temperature())
                   }
            domains_info.append(info)
        #__for domain in dna_structure.domain_list
        return domains_info 

    def _get_helices_info(self, dna_structure):
        """ Get JSON serialized data for helix objects. """
        helices_info = []
        # Need helices sorted by ID for indexing into helix arryay.
        helix_list = sorted(list(dna_structure.structure_helices_map.values()), key=lambda x: x.id)

        for helix in helix_list: 
            point1 = helix.end_coordinates[0]
            point2 = helix.end_coordinates[1]
            length = np.linalg.norm(point1-point2)
            frame = helix.end_frames[:,:,0]
            domain_ids = sorted(helix.get_domain_ids())
            connectivity_info = self._get_helix_conn_info(helix)
            possible_staple_crossovers = helix.possible_staple_crossovers
            possible_scaffold_crossovers = helix.possible_scaffold_crossovers

            info = { 'id'                 : helix.id,
                     'length'             : length,
                     'strand_radius'      : DnaParameters.strand_radius,
                     'base_pair_rise'     : DnaParameters.base_pair_rise,
                     'start_position'     : list(point1),
                     'orientation'        : [frame[0,2], frame[1,2], frame[2,2]],
                     'end_position'       : list(point2),
                     'scaffold_polarity'  : helix.scaffold_polarity,
                     'cadnano_info'       : { 'row' : helix.lattice_row, 'col' : helix.lattice_col, 
                                              'num' : helix.lattice_num },
                     'domains'            : domain_ids,
                     'helix_connectivity' : connectivity_info,
                     'num_possible_staple_crossovers' : len(possible_staple_crossovers),
                     'num_possible_scaffold_crossovers' : len(possible_scaffold_crossovers)
                   }

            helices_info.append(info)

        return helices_info

    def _get_strand_info(self, dna_structure):
        """ Get JSON serialized data for strands objetcs. 
        """
        strand_info_list = []
        for strand in dna_structure.strands:
            domain_ids = [domain.id for domain in strand.domain_list]
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
                     'domains'         : domain_ids,
                     'bases'           : base_info,
                     'color'           : strand.color
                   }

            strand_info_list.append(info)
        return strand_info_list

    def _get_helix_conn_info(self, helix):
        """ Get the information to write for helix connectivity.
        """
        #self._logger.setLevel(logging.DEBUG)
        self._logger.debug("==================== get conn information for helix num %d ===================" % helix.lattice_num) 
        dna_structure = self.dna_structure
        lattice = dna_structure.lattice
        num_neigh = lattice.number_of_neighbors
        self._logger.debug("Number of lattice directions: %d " % num_neigh) 
        helix_connectivity = helix.helix_connectivity 
        possible_staple_crossovers = helix.possible_staple_crossovers 
        possible_scaffold_crossovers = helix.possible_scaffold_crossovers
        self._logger.debug("Number of possible staple crossovers: %d " % len(possible_staple_crossovers)) 
        self._logger.debug("Number of possible scaffold crossovers: %d " % len(possible_scaffold_crossovers)) 
        self._logger.debug("Number of connected helices: %d " % len(helix_connectivity)) 

        # Create an array with size the number of possible helix neighors (4 for square lattice, 3 for honeycomb).
        helix_conn_array = [None]*num_neigh
        ang_inc = pi / num_neigh

        # Iterate over the helix connections filling in the appropriate element in helix_conn_array[]. 
        for connection in helix_connectivity: 
            to_helix = connection.to_helix
            nindex = lattice.get_neighbor_index(helix.lattice_row, helix.lattice_col, to_helix.lattice_row, to_helix.lattice_col)
            angle = ang_inc * nindex 
            dir = connection.direction
            num_staple = 0 
            num_scaffold = 0 
            crossovers = connection.crossovers
            for crossover in crossovers:
                if crossover.crossover_base.is_scaf:
                    num_scaffold += 1 
                else:
                    num_staple += 1 
            self._logger.debug(">>> Connected helix: num: %d  row: %d  col:%d " % (to_helix.lattice_num, to_helix.lattice_row,
                to_helix.lattice_col)) 
            self._logger.debug("    Nindex: %d " % nindex) 
            self._logger.debug("    Direction: (%g %g %g) " % (dir[0], dir[1], dir[2])) 

            crossover_info = self._get_crossover_info(connection)

            conn_info = { 'helix_id'   : to_helix.id,
                          'helix_num'  : to_helix.lattice_num,
                          'angle'      : angle,
                          'direction'  : list(dir),
                          'crossovers' : crossover_info
                        } 

            helix_conn_array[nindex] = conn_info 

        #__for connection in helix_connectivity

        #self._logger.debug(" Helix conn array: %s" % str(helix_conn_array)) 
        return helix_conn_array 

    def _get_crossover_info(self, connection):
        """ Get the helix crossover information to be written to a file. 
        """
        from_helix = connection.from_helix
        to_helix = connection.to_helix
        crossovers = connection.crossovers
        start_pos = from_helix.get_start_pos()
        self._logger.debug("    Start position: %d " % start_pos)

        # Create a list of staple and scaffold crossovers.
        staples = {} 
        scaffolds = {} 
        for crossover in crossovers:
            base = crossover.crossover_base
            if base.is_scaf:
                 scaffolds[base.p] = crossover
            else:
                staples[base.p] = crossover
        self._logger.debug("    Number of design crossovers: %d  staple: %d  scaffold: %d " % (len(crossovers), 
             len(staples), len(scaffolds))) 

        # Add staple information. Double crossovers will occur in pairs separated by 
        # a single position.
        crossover_info = [] 
        sorted_staples = collections.OrderedDict(sorted(staples.items()))
        self._logger.debug("    Staple crossovers: ") 
        pos = sorted_staples.keys()
        self._get_crossover_strand_info(start_pos, pos, staples, crossover_info)
        #print("########## staple crossover_info: size: %d" % len(crossover_info)) 
        #print("           data: %s" % str(crossover_info)) 

        # Add scaffold information.
        sorted_scaffolds = collections.OrderedDict(sorted(scaffolds.items()))
        self._logger.debug("    Scaffold crossovers: ") 
        pos = sorted_scaffolds.keys()
        self._get_crossover_strand_info(start_pos, pos, scaffolds, crossover_info)
        #print("########## scaffold crossover_info: size: %d" % len(crossover_info)) 
        #print("           data: %s" % str(crossover_info)) 

        return crossover_info 

    def _get_crossover_strand_info(self, start_pos, pos, crossovers, crossover_info):
        num_pos = len(pos)
        n = 0
        while n < num_pos:
            pos1 = pos[n]
            crossover1 = crossovers[pos1]
            base1 = crossover1.crossover_base
            strand1_id = crossover1.strand.id
            strand1_index = crossover1.strand.get_base_index(base1)
            add_pair = False
            if n != num_pos-1:
                pos2 = pos[n+1]
                if pos2-pos1 == 1:
                    crossover2 = crossovers[pos2]
                    base2 = crossover2.crossover_base
                    strand2_id = crossover2.strand.id
                    strand2_index = crossover2.strand.get_base_index(base2)
                    self._logger.debug("        n: %d  pos1: %d  pos2: %d " % (n, pos1, pos2))
                    self._logger.debug("               strand1: %d  index: %d " % (strand1_id, strand1_index))
                    self._logger.debug("               strand2: %d  index: %d " % (strand2_id, strand2_index))
                    add_pair = True 
                    n += 2

            if not add_pair:
                self._logger.debug("        n: %d  pos1: %d  str: %d" % (n, pos1, strand1_id))
                self._logger.debug("               strand1: %d  index: %d " % (strand1_id, strand1_index))
                pos2 = None
                strand2_id = None
                strand2_index = None 
                n += 1

            info = { "vhelix_base_index"        : pos1 - start_pos,
                     "first_strand_ID"          : strand1_id,
                     "first_strand_base_index"  : strand1_index,
                     "second_strand_ID"         : strand2_id,
                     "second_strand_base_index" : strand2_index
                   }

            crossover_info.append(info)
        #__while n < num_pos

        return crossover_info

    def _setup_logging(self):
        """ Set up logging."""
        self._logger = logging.getLogger(__name__)
        #self._logger = logging.getLogger('viewer:writer')
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

