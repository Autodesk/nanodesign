#!/usr/bin/env python
""" 
This module is used to write CanDo .cndo files. 
"""
import json
import logging
import numpy as np

class CandoWriter(object):
    """ The CandoWriter class writes out a CanDo .cndo file. 
    """
    def __init__(self, dna_structure):
        self.dna_structure = dna_structure
        self._logging_level = logging.INFO
        self._setup_logging()

    def set_logging_level(self,level):
        """Set logging level."""
        self._logger.setLevel(level)

    def write(self,file_name):
        """Write a .cndo file.

        Args:
            file_name (string): The name of a viewer JSON file to write.

        """
        dna_structure = self.dna_structure
        base_conn = dna_structure.base_connectivity 
        base_nodes = dna_structure.helix_axis_nodes
        triad = dna_structure.helix_axis_frames
        id_nt = dna_structure.id_nt
        self._logger.info("Writing CanDo .cndo file: %s " % file_name)
        self._logger.info("Number of bases in the base connectivity table: %d " % len(base_conn))
        self._logger.info("Number of base helix axis nodes: %d " % len(base_nodes))

        with open(file_name, 'w') as cndo_file:
            # write header
            cndo_file.write('"CanDo (.cndo) file format version 1.0, Keyao Pan, Laboratory for Computational Biology and Biophysics, Massachusetts Institute of Technology, November 2015"\n');
            cndo_file.write("\n")

            # write dna topology
            cndo_file.write("dnaTop,id,up,down,across,seq\n")
            for i in xrange(0,len(base_conn)):
                base = base_conn[i] 
                cndo_file.write("%d,%d,%d,%d,%d,%s\n" % (i+1, base.id, base.up, base.down, base.across, base.seq))
            cndo_file.write("\n")

            # base nodes
            cndo_file.write('dNode,"e0(1)","e0(2)","e0(3)"\n')
            for i in xrange(0,len(base_nodes)):
                cndo_file.write("%d,%f,%f,%f\n" % (i+1, base_nodes[i,0], base_nodes[i,1], base_nodes[i,2]))
            cndo_file.write("\n")

            # triad vectors
            cndo_file.write('triad,"e1(1)","e1(2)","e1(3)","e2(1)","e2(2)","e2(3)","e3(1)","e3(2)","e3(3)"\n')
            for i in xrange(0,triad.shape[2]):
                cndo_file.write("%d,%f,%f,%f,%f,%f,%f,%f,%f,%f\n" % (i+1, -triad[0,0,i], -triad[1,0,i], -triad[2,0,i],
                    triad[0,1,i], triad[1,1,i], triad[2,1,i], -triad[0,2,i], -triad[1,2,i], -triad[2,2,i]))
            cndo_file.write("\n")

            # nucleotide binding table 
            cndo_file.write("id_nt,id1,id2\n")
            for i in xrange(0,len(id_nt)):
                cndo_file.write("%d,%d,%d\n" % (i+1, id_nt[i,0]+1, id_nt[i,1]+1))

        self._logger.info("Done.")

    def _setup_logging(self):
        """ Set up logging."""
        self._logger = logging.getLogger('cando:writer')
        self._logger.setLevel(self._logging_level)
        # create console handler and set format
        console_handler = logging.StreamHandler()
        formatter = logging.Formatter('[%(name)s] %(levelname)s - %(message)s')
        console_handler.setFormatter(formatter)
        self._logger.addHandler(console_handler)

def main():
    """ Write a CanDo .cndo file."""
    file_name = sys.argv[1]
    cadnano_reader = CadnanoReader()
    #cadnano_reader.set_logging_level(logging.DEBUG)
    cadnano_reader.read_json(file_name)

if __name__ == '__main__':
    main()

