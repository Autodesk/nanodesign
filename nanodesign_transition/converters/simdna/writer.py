#!/usr/bin/env python
""" 
This module is used to write SimDNA pairs files. 
"""
import collections
import logging

try:
    import os.path
    import sys
    base_path = os.path.abspath( os.path.dirname(__file__) + '/../../../' )
    sys.path.append( base_path )
    sys.path = sys.path[:-1]
except ImportError as i:
    print "Could not get nanodesign_transition module"
    raise i

class SimDnaWriter(object):
    """ The SimDnaWriter class writes out a SimDNA pairs file. 
    """
    def __init__(self, dna_structure):
        self.dna_structure = dna_structure
        self._logging_level = logging.INFO
        self._setup_logging()

    def set_logging_level(self,level):
        """Set logging level."""
        self._logger.setLevel(level)

    def write(self,file_name):
        """Write a pairs file.

        Args:
            file_name (string): The name of pairs file to write.

        """

        self._logger.info("Writing SimDNA pairs file %s " % file_name)
        dna_structure = self.dna_structure
        dna_structure.compute_aux_data()
        num_bases = len(dna_structure.base_connectivity)
        nm_to_ang = 10.0

        with open(file_name, 'w') as outfile:
            outfile.write("%d\n" % num_bases)
            for strand in dna_structure.strands:
                base_coords = strand.get_base_coords()
                for i in xrange(0,len(strand.tour)):
                    base = strand.tour[i]
                    if base.across == None or base.across == -1:
                        paired_strand_id = -1
                        paired_base_id = -1
                    else:
                        across_base = base.across
                        paired_strand_id = across_base.strand
                        paired_strand = self.dna_structure.strands_map[paired_strand_id]
                        #print(">>> strand id %d  paired_strand  id %d" % (strand.id, paired_strand.id))
                        paired_base_id = paired_strand.get_base_index(across_base)+1
                    #coord = nm_to_ang * base_coords[i]
                    coord = nm_to_ang * base.coord
                    base_id = strand.get_base_index(base)+1
                    outfile.write("%4d %4d %8g %8g %8g %4d %4d\n" % 
                        (strand.id, base.id, coord[0], coord[1], coord[2], paired_strand_id, paired_base_id))
                #__for i in xrange(0,len(strand.tour))
            #__for strand in dna_structure.strands
        #__with open(file_name, 'w') as outfile

    def _setup_logging(self):
        """ Set up logging."""
        self._logger = logging.getLogger(__name__)
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

