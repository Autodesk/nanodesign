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
        self._logger = logging.getLogger(__name__)   

    def write(self,file_name):
        """Write a .cndo file.

        Args:
            file_name (string): The name of a viewer JSON file to write.

        """
        dna_structure = self.dna_structure
        base_conn = dna_structure.base_connectivity 
        self._logger.info("Writing CanDo .cndo file: %s " % file_name)
        self._logger.info("Number of bases %d " % len(base_conn))

        with open(file_name, 'w') as cndo_file:
            # write header
            cndo_file.write('"CanDo (.cndo) file format version 1.0, Keyao Pan, Laboratory for Computational Biology and Biophysics, Massachusetts Institute of Technology, November 2015"\n');
            cndo_file.write("\n")

            # write dna topology
            cndo_file.write("dnaTop,id,up,down,across,seq\n")
            for i in xrange(0,len(base_conn)):
                base = base_conn[i] 
                up = base.up.id if base.up else -1
                down = base.down.id if base.down else -1
                across = base.across.id if base.across else -1
                cndo_file.write("%d,%d,%d,%d,%d,%s\n" % (i+1, base.id, up, down, across, base.seq))
            cndo_file.write("\n")

            # base nodes
            cndo_file.write('dNode,"e0(1)","e0(2)","e0(3)"\n')
            for i in xrange(0,len(base_conn)):
                coords = base_conn[i].coordinates 
                cndo_file.write("%d,%f,%f,%f\n" % (i+1, coords[0], coords[1], coords[2]))
            cndo_file.write("\n")

            # triad vectors
            cndo_file.write('triad,"e1(1)","e1(2)","e1(3)","e2(1)","e2(2)","e2(3)","e3(1)","e3(2)","e3(3)"\n')
            for i in xrange(0,len(base_conn)):
                ref_frame = base_conn[i].ref_frame
                cndo_file.write("%d,%f,%f,%f,%f,%f,%f,%f,%f,%f\n" % (i+1, -ref_frame[0,0], -ref_frame[1,0], -ref_frame[2,0],
                    ref_frame[0,1], ref_frame[1,1], ref_frame[2,1], -ref_frame[0,2], -ref_frame[1,2], -ref_frame[2,2]))
            cndo_file.write("\n")

            # Nucleotide binding table.
            id_nt = self._create_id_nt(base_conn)
            cndo_file.write("id_nt,id1,id2\n")
            for i in xrange(0,len(id_nt)):
                cndo_file.write("%d,%d,%d\n" % (i+1, id_nt[i][0]+1, id_nt[i][1]+1))
        self._logger.info("Done.")

    def _setup_logging(self):
        """ Set up logging."""
        self._logger = logging.getLogger(__name__)
        self._logger.setLevel(self._logging_level)
        # Create console handler and set format.
        if not len(self._logger.handlers):
            console_handler = logging.StreamHandler()
            formatter = logging.Formatter('[%(name)s] %(levelname)s - %(message)s')
            console_handler.setFormatter(formatter)
            self._logger.addHandler(console_handler)

    def _create_id_nt(self, base_conn):
        """ Create a list of paired bases. """
        id_nt = []
        for i,base in enumerate(base_conn):
            if not (base.across and base.is_scaf):
                continue
            id_nt.append([base.id, base.across.id])
        #__for i,base in enumerate(base_conn):
        return id_nt
    #__def _create_nt_id

