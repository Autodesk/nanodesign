#!/usr/bin/env python
""" 
This module is used to read PDB format files used internally to create atomic models of DNA nanodesigns. 

The PDB files are used as templates to create atomic models for nanodesigns. Only ATOM records are read.
"""
import logging
import sys

import atomic_structure
# TODO (JS 10/26/16): There's a weird circular dependency here still. We got
# around it due to importing the names at the same level, but atomic_structure
# needs PdbReader and PdbReader needs Atom and Molecule.

# PDB ATOM record format. 
ATOM_FORMAT = "%6s%5d %4s%c%3s %c%4d%c   %8.3f%8.3f%8.3f%6.2f%6.2f           %2s \n"

class PdbRecordTypes(object):
    """ PDB record names. """
    ATOM = "ATOM"
    TER  = "TER"
    names   = [ ATOM, TER ]

class PdbReader(object):
    """ The PdbReader class reads in a PDB format .pdb file. 
    """
    logging_set = False

    def __init__(self):
        #self._logging_level = logging.INFO
        #self._setup_logging()
        self.current_molecule = None
        self.molecules = []
        self.process_record = { PdbRecordTypes.ATOM : PdbReader.process_atom,
                                PdbRecordTypes.TER  : PdbReader.process_ter
                              }
        self._logging_level = logging.INFO
        self._setup_logging()

    def set_logging_level(self,level):
        """Set logging level."""
        self._logger.setLevel(level)

    def _setup_logging(self):
        """ Set up logging."""
        self._logger = logging.getLogger(__name__)
        self._logger.setLevel(self._logging_level)
        if PdbReader.logging_set:
           return
        # Create console handler and set format.
        if not len(self._logger.handlers):
            console_handler = logging.StreamHandler()
            formatter = logging.Formatter('[%(name)s] %(levelname)s - %(message)s')
            console_handler.setFormatter(formatter)
            self._logger.addHandler(console_handler)
            PdbReader.logging_set = True

    def read(self,file_name):
        """Read a .pdb file.

        Args:
            file_name (string): The name of the PDB file to read.
        """
        self._logger.info("Reading PDB file %s " % file_name)

        # Create an initial molecule to read data in to.
        molecule = atomic_structure.Molecule(1)
        self.current_molecule = molecule
        self.molecules.append(molecule)

        # Read the file processing the PDB records.
        with open(file_name, "r") as pdb_file:
            for line in pdb_file:
               rec_type = line[0:4].strip()
               if (rec_type in self.process_record):
                   self.process_record[rec_type](self,line)
               else:
                   self._logger.warn("Unknown record type: %s " % rec_type)
        #_with open_

        self._logger.info("Done.")
        self._logger.info("Read %d atoms." % len(self.current_molecule.atoms))
        self._logger.info("Read %d chains: %s " % (len(self.current_molecule.chains), str(list(self.current_molecule.chains))))

    def process_atom(self,line):
        """ Process a PDB ATOM record. """
        atom_id = int(line[6:11])
        atom_name = line[12:16].strip()
        alt_loc = line[16]
        res_name = line[17:20].strip()
        chain_id = line[21]
        seq_id = int(line[22:26])
        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])
        element = line[76:78]
        atom = atomic_structure.Atom(atom_id, atom_name, res_name, chain_id, seq_id, x, y, z, element )
        self.current_molecule.add_atom(atom)
    #__def process_atom

    def process_ter(self,line):
        """ Process a PDB TER record. """
        pass
    #__def process_ter

