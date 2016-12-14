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
This module is used to write CIF formated files for the atomic structure of a dna structure. 

The atomic structure is generated using an AtomicStructure object and an input DnaStructure object.
A Molecular object is created for each strand in the dna structure and contains the atoms for each
base residue the strand.

A single entity is created (id=1) with the number of chains equal to the number of strands in the dna structure.
The chain IDs for each strand are generated using the format: 

    <sc|st>.<vhelixNum>.<startPos>

    where sc=scaffold, st=staple
          vhelixNUm = the virtual helix number from cadnano
          startPos = the position in the virtual helix of the first base in the strand.

"""
from collections import OrderedDict
import json
import logging
import os 
import numpy as np
from .atomic_structure import AtomicStructure

class CifWriter(object):
    """ The CifWriter class is used to write a CIF formated file. 

        Attributes:
            dna_structure (DnaStructure) : The dna structure to convert to an atomic structure and write to a CIF file.
            entityID (int): The CIF entiy ID for the dna structure.
    """

    # Define some constants used in the CIF file.
    STRUCT_ENTRY_ID = "Nanodesign"
    ENTITY_PDBX_DESCRIPTION = '"Nanodesign structure"'
    EMPTY_FIELD = "?"
    COMMENT_SPACES = "#\n#\n"
    LOOP = "loop_\n"

    def __init__(self, dna_structure):
        """
            Initialize the CifWriter object.

            Arguments:
                dna_structure (DnaStructure) : The dna structure to convert to an atomic structure and write to a CIF file.
        """
        self.dna_structure = dna_structure
        self.entityID = 1
        self._logger = logging.getLogger(__name__)   

    def write(self, file_name, infile, informat):
        """ Write a CIF file.

            Arguments:
                file_name (string): The name of the CIF file to write.
                infile (string): The name of the file the DNA structure was created from.
                informat (string): The format of the file the DNA structure was created from.
        """
        dna_structure = self.dna_structure
        base_conn = dna_structure.base_connectivity
        strands = dna_structure.strands
        self._logger.info("Writing CIF file %s " % file_name)
        self._logger.info("Number of bases %d " % len(base_conn))
        self._logger.info("Number of strands %d " % len(strands))

        # Generate atomic models of the dna structure. A list of Molecule objects is 
        # created for each strand.  
        atomic_structure = AtomicStructure(dna_structure)
        molecules = atomic_structure.generate_structure_ss()  # converts ssDNA
        #molecules = atomic_structure.generate_structure()
        self._logger.info("Number of molecules %d " % len(molecules))
        num_atoms = 0
        for molecule in molecules:
            num_atoms += len(molecule.atoms)
        self._logger.info("Number of atoms %d " % num_atoms) 

        # Write the CIF file.
        with open(file_name, 'w') as cif_file:
            cif_file.write("data_nanodesign_structure")
            cif_file.write(CifWriter.COMMENT_SPACES)
            self._write_struct_records(cif_file, informat, infile)
            self._write_entity_records(cif_file, atomic_structure.strands)
            self._write_struct_asym_records(cif_file, atomic_structure.strands)
            self._write_atom_site_records(cif_file, molecules)
        self._logger.info("Done.")

    def _write_struct_records(self, cif_file, informat, infile):
        """ Write CIF _struct records. 

            Arguments:
                cif_file (file): The file object used to write the CIF file.
                infile (string): The name of the file the DNA structure was created from.
                informat (string): The format of the file the DNA structure was created from.
        """
        _,infile_name = os.path.split(infile)

        # Define an OrderedDict of field names and values for the _struct records. 
        struct_records = OrderedDict([
            ("entry_id"                , CifWriter.STRUCT_ENTRY_ID), 
            ("title"                   , '"CIF file generated from %s format file %s"' % (informat, infile_name)), 
            ("pdbx_descriptor"         , CifWriter.EMPTY_FIELD), 
            ("pdbx_model_details"      , CifWriter.EMPTY_FIELD),
            ("pdbx_CASP_flag"          , CifWriter.EMPTY_FIELD),
            ("pdbx_model_type_details" , CifWriter.EMPTY_FIELD) 
        ])

        for name,value in struct_records.iteritems(): 
            cif_file.write("_struct.%-30s       %s\n" % (name, value)) 
        cif_file.write(CifWriter.COMMENT_SPACES)

    def _write_entity_records(self, cif_file, strands):
        """ Write CIF _entity records. 

            Arguments:
                cif_file (file): The file object used to write the CIF file.
                strands (List[AtomicStructureStrand]): The list of the strands for the DNA atomic structure. 

            A single entity is defined with ID 1. The number of molecules is equal to the number of strands in the DNA structure.
        """

        # Define an OrderedDict of field names and values for the _entity records. 
        entity_records = OrderedDict([
            ("id"                      , self.entityID),
            ("type"                    , "polymer"),
            ("pdbx_description"        , CifWriter.ENTITY_PDBX_DESCRIPTION),   
            ("formula_weight"          , CifWriter.EMPTY_FIELD),
            ("pdbx_number_of_molecule" , len(strands)),
            ("details", CifWriter.EMPTY_FIELD)
        ])

        # Write the _entity records.
        for name,value in entity_records.iteritems():
            cif_file.write("_entity.%-30s       %s\n" % (name, value))
        cif_file.write(CifWriter.COMMENT_SPACES)

    def _write_struct_asym_records(self, cif_file, strands):
        """ Write CIF _struct_asym records. 

            Arguments:
                cif_file (file): The file object used to write the CIF file.
                strands (List[AtomicStructureStrand]): The list of the strands for the DNA atomic structure. The chain IDs for the
                   _struct_asym records are obtained from the strands in this list.
          
            The _struct_asym records define the chain IDs for each of the entities in the structure.
        """

        # Define the list of field names for the _struct_asym records. 
        struct_asym_records = [ "id", "pdbx_blank_PDB_chainid_flag", "pdbx_modified", "entity_id", "details" ]
        cif_file.write(CifWriter.LOOP)

        # Write the _struct_asym record field names.
        for name in struct_asym_records:
            cif_file.write("_struct_asym.%s\n" % name)

        # Write the _struct_asym record data.
        for strand in strands:
            cif_file.write("%s N N %d %s\n" % (strand.chainID, self.entityID, CifWriter.EMPTY_FIELD))
        cif_file.write(CifWriter.COMMENT_SPACES)

    def _write_atom_site_records(self, cif_file, molecules):
        """ Write CIF _atom_site records. 

            Arguments:
                cif_file (file): The file object used to write the CIF file.
                molecules(List[AtomicStructureMolecule]): The list of the molecules for the DNA atomic structure. 

            The CIF _atom_site records describe the atom data for the DNA structure. The atom data is obtained 
            from a Molecule object.
        """

        # Define an OrderedDict of field names,values and print format. 
        # Note (davep) I was using this for a nice (but slow) way to keep track of data
        # and fields. I'll keep it around for now because it gives a nice visual for field data.
        atom_site_records = OrderedDict([
            ("group_PDB",         ["ATOM", "%4s"]),
            ("id",                [None, "%6d"]),
            ("type_symbol",       [None, "%5s"]),
            ("label_atom_id",     [None, "%3s"]),
            ("label_alt_id",      [".", "%4s"]),
            ("label_comp_id",     [None, "%s"]),
            ("label_asym_id",     [None, "%s"]),
            ("label_entity_id",   [self.entityID, "%d"]),
            ("label_seq_id",      [None, "%5d"]),
            ("pdbx_PDB_ins_code", [CifWriter.EMPTY_FIELD, "%s"]),
            ("Cartn_x",           [None, "%g"]),
            ("Cartn_y",           [None, "%g"]),
            ("Cartn_z",           [None, "%g"]),
            ("occupancy",         [1.0, "%g"]),
            ("B_iso_or_equiv",    [1.0, "%g"]),
            ("Cartn_x_esd",       [CifWriter.EMPTY_FIELD, "%s"]),
            ("Cartn_y_esd",       [CifWriter.EMPTY_FIELD, "%s"]),
            ("Cartn_z_esd",       [CifWriter.EMPTY_FIELD, "%s"]),
            ("occupancy_esd",     [CifWriter.EMPTY_FIELD, "%s"]),
            ("B_iso_or_equiv_esd",[CifWriter.EMPTY_FIELD, "%s"]),
            ("pdbx_formal_charge",[CifWriter.EMPTY_FIELD, "%s"]),
            ("auth_seq_id",       [None, "%d"]),
            ("auth_comp_id",      [None, "%s"]),
            ("auth_asym_id",      [None, "%s"]),
            ("auth_atom_id",      [None, "%s"]),
            ("pdbx_PDB_model_num",[1,  "%d"])
       ])
  
        # Set the values of fields that doe not change.
        label_alt_id = "."
        label_entity_id = self.entityID
        pdbx_PDB_ins_code = CifWriter.EMPTY_FIELD
        occupancy = 1.0
        B_iso_or_equiv = 1.0
        Cartn_x_esd = CifWriter.EMPTY_FIELD
        Cartn_y_esd = CifWriter.EMPTY_FIELD
        Cartn_z_esd = CifWriter.EMPTY_FIELD
        occupancy_esd = CifWriter.EMPTY_FIELD
        B_iso_or_equiv_esd = CifWriter.EMPTY_FIELD
        pdbx_formal_charge = CifWriter.EMPTY_FIELD
        pdbx_PDB_model_num = 1

        # Create the atom record format.
        format = ""
        for _,value in atom_site_records.iteritems():
            format += value[1] + "  "
        format += "\n"

        # Write atom field names.
        cif_file.write(CifWriter.LOOP)
        for name,value in atom_site_records.iteritems():
            cif_file.write("_atom_site.%s\n" % name)

        # Write atom data.
        for molecule in molecules:
            for atom in molecule.atoms:
                id = atom.id
                type_symbol = atom.element
                label_atom_id = atom.name
                label_comp_id = atom.res_name
                label_asym_id = atom.chainID
                label_seq_id = atom.res_seq_num
                Cartn_x = atom.coords[0]
                Cartn_y = atom.coords[1]
                Cartn_z = atom.coords[2]
                auth_seq_id = atom.res_seq_num
                auth_comp_id = atom.res_name
                auth_asym_id = atom.chainID
                auth_atom_id = atom.name
                cif_file.write(format % ("ATOM ", id, type_symbol, label_atom_id, label_alt_id, label_comp_id, label_asym_id, 
                    label_entity_id, label_seq_id, pdbx_PDB_ins_code, Cartn_x, Cartn_y, Cartn_z, occupancy, B_iso_or_equiv, 
                    Cartn_x_esd, Cartn_y_esd, Cartn_z_esd, occupancy_esd, B_iso_or_equiv_esd, pdbx_formal_charge, auth_seq_id, 
                    auth_comp_id, auth_asym_id, auth_atom_id, pdbx_PDB_model_num))
            #__for atom in molecule.atoms
        #__for molecule in molecules


