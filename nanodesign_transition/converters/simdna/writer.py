#!/usr/bin/env python
""" 
This module is used to write SimDNA pairs files. 

The SimDNA pairs file file format is: 

    1st line is just an integer, indicating number of base records in the file.

    Remaining lines are 7-columns space separated records: 

        1st column: Chain (strand) id. Integer, 0-indexed

        2nd column: Base id. Integer, 1-indexed. Strand-relative

        3rd-5th column: X,Y,Z coordinates of C5' atom, floating point.

        6th column: Chain (strand) id of paired base. Integer, 0-indexed. -1 if base is unpaired

        7th column: Base id of paired base. Integer, 1-indexed. Strand-relative. -1 if base is unpaired.

    The base records are sorted by chain (ascending), and then by base id within that chain (ascending).
    The chain with id 0 is the scaffold. It's unclear how this will work with multiple scaffolds, but 0 
    should always be the scaffold if possible.

    All bases are indexed using strand-relative base ids. So a strand with length 32 will have base ids 1-32, etc.
"""
import collections
import itertools
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
        """ Write a pairs file.

        Arguments:
            file_name (string): The name of pairs file to write.
        """
        self._logger.info("Writing SimDNA pairs file %s " % file_name)
        dna_structure = self.dna_structure
        dna_structure.compute_aux_data()
        num_bases = len(dna_structure.base_connectivity)
        nm_to_ang = 10.0

        with open(file_name, 'w') as outfile:
            outfile.write("%d\n" % num_bases)
            # Create a list of scaffold and staple strands.
            scaffold_strands = []
            staple_strands = []
            for strand in dna_structure.strands:
                if strand.is_scaffold:
                    scaffold_strands.append(strand)
                else:
                    staple_strands.append(strand)
            #__for strand in dna_structure.strands

            # Redefine output strand IDs to make sure scaffold strand has a 0 ID.
            strand_id = 0
            strand_map = {}
            for strand in itertools.chain(scaffold_strands, staple_strands):
                strand_map[strand.id] = strand_id
                strand_id += 1
            #__for strand in itertools.chain(scaffold_strands, staple_strands)

            # Write base records.
            for strand in itertools.chain(scaffold_strands, staple_strands):
                strand_id = strand_map[strand.id]
                base_coords = strand.get_base_coords()

                for i in xrange(0,len(strand.tour)):
                    base = strand.tour[i]
                    if base.across == None:
                        paired_strand_id = -1
                        paired_base_id = -1
                    else:
                        across_base = base.across
                        paired_strand_id = across_base.strand
                        paired_strand = self.dna_structure.strands_map[paired_strand_id]
                        paired_base_id = paired_strand.get_base_index(across_base)+1
                    #__if base.across == None
                    coord = nm_to_ang * base.nt_coords
                    base_id = strand.get_base_index(base)+1
                    outfile.write("%4d %4d %8g %8g %8g %4d %4d\n" % 
                        (strand_id, i+1, coord[0], coord[1], coord[2], paired_strand_id, paired_base_id))
                #__for i in xrange(0,len(strand.tour))
            #__for strand in itertools.chain(scaffold_strands, staple_strands)
        #__with open(file_name, 'w') as outfile

    def _setup_logging(self):
        """ Set up logging."""
        self._logger = logging.getLogger(__name__)
        self._logger.setLevel(self._logging_level)

        # create console handler and set format
        if not len(self._logger.handlers):
            console_handler = logging.StreamHandler()
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

