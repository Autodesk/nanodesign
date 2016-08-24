#!/usr/bin/en
""" 
This module is used to write atomic structures to a file in the PDB format. 

Only the ATOM records are written. Each strand is stored within a PDB MODEL record. This allows the
PDB file to be visualized using Chimera. Strands are also given a single letter or digit chain ID take n
from [A-Z,a-z,0-9], cycling through the list as needed. 
"""
import json
import logging
import numpy as np
from math import sqrt 
from nanodesign_transition.atomic_structure import AtomicStructure 

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

class PdbWriter(object):
    """ The PdbWriter class writes out a PDB format .pdb file. 
    """
    # Formats used to write PDB records.
    ATOM_FORMAT = "ATOM %6d %4s%1s%3s%2s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f      %4s%2s%2s\n"
    MODEL_FORMAT = "MODEL     %4d\n"
    ENDMDL_FORMAT = "ENDMDL\n"
    TER_FORMAT = "TER  %6d      %3s%2s%4d%1s\n" 

    # A list of single letter or digit chain IDs used for strands.
    CHAIN_IDS = list('ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789')

    def __init__(self, dna_structure):
        self.dna_structure = dna_structure
        self._logging_level = logging.INFO
        self._setup_logging()

    def set_logging_level(self,level):
        """Set logging level."""
        self._logger.setLevel(level)

    def write(self,file_name):
        """Write a .pdb file.

        Arguments:
            file_name (string): The name of the PDB file to write.
        """
        dna_structure = self.dna_structure
        base_conn = dna_structure.base_connectivity
        strands = dna_structure.strands
        triad = dna_structure.helix_axis_frames
        id_nt = dna_structure.id_nt
        self._logger.info("Writing PDB .pdb file: %s " % file_name)
        self._logger.info("Number of bases %d " % len(base_conn))
        self._logger.info("Number of strands %d " % len(strands))

        # Generate atomic models of the dna structure.
        atomic_structure = AtomicStructure(dna_structure)
        molecules = atomic_structure.generate_structure_ss()   # converts ssDNA 
        #molecules = atomic_structure.generate_structure()
        xmin,xmax,ymin,ymax,zmin,zmax = atomic_structure.get_extent()

        # Write the models.
        self.set_logging_level(logging.DEBUG)
        res_seq = 0 
        model_num = 1
        atom_id = 1
        chain_count = 0
        with open(file_name, 'w') as pdb_file:
            # write header
            #pdb_file.write('"PDB file generated from "');
            pdb_file.write(PdbWriter.MODEL_FORMAT % model_num)
            for molecule in molecules:
                chain_id = PdbWriter.CHAIN_IDS[chain_count] 
                res_seq, atom_id, model_num = self._write_molecule(pdb_file, molecule, model_num, res_seq, atom_id, 
                    chain_id, xmin, ymin, zmin)
                #model_num += 1
                chain_count += 1
                if chain_count == len(PdbWriter.CHAIN_IDS):
                    chain_count = 0

            pdb_file.write("\n")
            pdb_file.write(PdbWriter.ENDMDL_FORMAT)
        self._logger.info("Done.")

    def _setup_logging(self):
        """ Set up logging."""
        self._logger = logging.getLogger(__name__)
        self._logger.setLevel(self._logging_level)
        # create console handler and set format
        console_handler = logging.StreamHandler()
        formatter = logging.Formatter('[%(name)s] %(levelname)s - %(message)s')
        console_handler.setFormatter(formatter)
        self._logger.addHandler(console_handler)

    def _write_molecule(self, pdb_file, molecule, model_num, res_seq, atom_id, chain_id, xmin, ymin, zmin):
        """ Write the atoms in a molecule to a file. 

            Arguments:
                pdb_file (File): The file handle used to write to the file. 
                molecule (Molecule): The Molecule object to write.
                model_num (int): 
                res_seq (int): 
                atom_id (int):, 
                chain_id (string): 
                xmin (Float): The x minimum extent of the atomic structure.
                ymin (Float): The y minimum extent of the atomic structure.
                zmin (Float): The z minimum extent of the atomic structure.
        """
        self._logger.debug("Write molecule %d " % molecule.model_id)
        self._logger.debug("Number of chains %d  %s" % (len(molecule.chains), str(list(molecule.chains))))
        self._logger.debug("Number of atoms %d " % (len(molecule.atoms)))
        cmpd = "   "
        current_res_seq = molecule.atoms[0].res_seq_num
        res_seq += 1 
        #pdb_file.write(PdbWriter.MODEL_FORMAT % model_num)
        new_res = False

        for atom in molecule.atoms:
            x = atom.coords[0] - xmin
            y = atom.coords[1] - ymin
            z = atom.coords[2] - zmin
            #id = atom.id
            id = atom_id
            name = atom.name.ljust(3)
            res = atom.res_name
            if current_res_seq != atom.res_seq_num:
                current_res_seq = atom.res_seq_num
                res_seq += 1
                new_res = True
            else:
                new_res = False

            seq = atom.res_seq_num
            chain = chain_id
            element = atom.element
            mass = 1.0
            remote_ind = " "
            branch = " "
            alt_loc = " "
            icode = " "
            occupancy = 0.0
            temp_factor = 0.0
            segID = " "
            charge = " "

            if atom_id > 99900 and new_res: 
                pdb_file.write(PdbWriter.ENDMDL_FORMAT)
                model_num += 1
                atom_id = 1 
                id = atom_id 
                pdb_file.write(PdbWriter.MODEL_FORMAT % model_num)

            pdb_file.write(PdbWriter.ATOM_FORMAT % (id, name, alt_loc, res, chain, seq, icode, x, y, z, occupancy, 
                temp_factor, segID, element, charge))
            atom_id += 1

        #__for atom in molecule.atoms

        pdb_file.write(PdbWriter.TER_FORMAT % (atom_id, res, chain, seq, icode))
        atom_id += 1
        #pdb_file.write(PdbWriter.ENDMDL_FORMAT)
        return res_seq, atom_id, model_num

