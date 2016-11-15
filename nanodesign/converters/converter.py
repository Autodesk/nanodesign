#!/usr/bin/env python
""" 
This module is used to convert DNA design files into other file formats.

An input caDNAno design file is conveted into a DnaStructure object containing information about the
design (e.g. lattice type, virtual helix definitions) and information derived from that design 
(e.g. strands, domains). caDNAno design files may contain deleted/inserted bases. By default the 
DnaStructure is not created with deleted/inserted bases. The DnaStructure is created with 
deleted/inserted bases by specifying the --modify command-line argument.

"""
import os
import re
import sys
import json
import logging
from .cadnano.reader import CadnanoReader
from .cadnano.writer import CadnanoWriter
from .cadnano.convert_design import CadnanoConvertDesign
from .viewer.writer import ViewerWriter 
from .cando.writer import CandoWriter 
from .simdna.writer import SimDnaWriter 
from .pdbcif.pdb_writer import PdbWriter 
from .pdbcif.cif_writer import CifWriter 

from ..data.dna_structure import DnaStructure
from ..data.parameters import DnaParameters
from ..utils.xform import Xform,HelixGroupXform,apply_helix_xforms,xform_from_connectors

from .dna_sequence_data import dna_sequence_data
# TODO (JMS, 10/26/16): revisit where the sequence data is kept?

class ConverterFileFormats(object):
    """ File format names to convert to/from. """
    UNKNOWN   = "unknown"
    CADNANO   = "cadnano"
    CANDO     = "cando"
    CIF       = "cif"
    PDB       = "pdb"
    SIMDNA    = "simdna"
    STRUCTURE = "structure"
    TOPOLOGY  = "topology"
    VIEWER    = "viewer"
    names = [ CADNANO, CANDO, CIF, PDB, SIMDNA, STRUCTURE, TOPOLOGY, VIEWER ]

class Converter(object):
    """ This class stores objects for various models created when reading from a file.

        Attributes:
            cadnano_design (CadnanoDesign): The object storing the caDNAno design information.
            cadnano_convert_design (CadnanoConvertDesign): The object used to convert a caDNAno design into a DnaStructure.
            dna_parameters (DnaParameters): The DNA physical parameters used to generate the geometry of a DNA structure
            dna_structure (DnaStructure): The object storing connectivity and geometry of a DNA structure. 
            infile (String): The file name to convert.
            informat (String): The format of the file to convert, taken from ConverterFileFormats.
            modify (bool): If true then DnaStructure is created with deleted/inserted bases.
            outfile (String): The name of the file for converter output.
    """
    def __init__(self):
        self.cadnano_design = None 
        self.dna_structure = None 
        self.cadnano_convert_design = None
        self.infile = None
        self.informat = None
        self.outfile = None
        self.modify = False
        self.dna_parameters = DnaParameters()
        self.logger = logging.getLogger(__name__)

    def read_cadnano_file(self, file_name, seq_file_name, seq_name):
        """ Read in a caDNAno file. 

            Arguments:
                file_name (String): The name of the caDNAno file to convert. 
                seq_file_name (String): The name of the CSV file used to assign a DNA base sequence to the DNA structure. 
                seq_name (String): The name of a sequence used to assign a DNA base sequence to the DNA structure.
        """
        cadnano_reader = CadnanoReader()
        self.cadnano_design = cadnano_reader.read_json(file_name)
        self.cadnano_convert_design = CadnanoConvertDesign(self.dna_parameters)
        self.dna_structure = self.cadnano_convert_design.create_structure(self.cadnano_design, self.modify)

        # Read in staple sequences from a CSV format file.
        if (seq_file_name): 
            _, file_extension = os.path.splitext(seq_file_name)

            if (file_extension == ".csv"): 
                modified_structure = False
                sequence = cadnano_reader.read_csv(seq_file_name)
                self.cadnano_convert_design.set_sequence(self.dna_structure, modified_structure, sequence)

        # Assign a sequence using a name.
        if (seq_name): 
            if (seq_name not in dna_sequence_data):
                self.logger.error("The sequence name %s is not recognized.", seq_name)
            modified_structure = False
            self.cadnano_convert_design.set_sequence_from_name(self.dna_structure, modified_structure, seq_name)

    def write_viewer_file(self, file_name):
        """ Write a Nanodesign Viewer file.

            Arguments:
                file_name (String): The name of the Nanodesign Viewer file to write. 
        """
        viewer_writer = ViewerWriter(self.dna_structure, self.dna_parameters)
        viewer_writer.write(file_name)

    def write_pdb_file(self, file_name):
        """ Write an RCSB PDB-format file.

            Arguments:
                file_name (String): The name of the PDB file to write. 
        """
        pdb_writer = PdbWriter(self.dna_structure)
        pdb_writer.write(file_name)

    def write_cif_file(self, file_name):
        """ Write a RCSB CIF-format file.

            Arguments:
                file_name (String): The name of the CIF file to write. 
        """
        cif_writer = CifWriter(self.dna_structure)
        cif_writer.write(file_name, self.infile, self.informat )

    def write_simdna_file(self, file_name):
        """ Write a SimDNA pairs file.

            Arguments:
                file_name (String): The name of the SimDNA pairs file to write. 
        """
        simdna_writer = SimDnaWriter(self.dna_structure)
        simdna_writer.write(file_name)

    def write_topology_file(self, file_name):
        """ Write a DNA topology file.

            Arguments:
                file_name (String): The name of the topology file to write. 
        """
        self.dna_structure.write_topology(file_name, write_json_format=True)

    def write_structure_file(self, file_name):
        """ Write a DNA structure file.
            Arguments:
                file_name (String): The name of the structure file to write. 
        """
        self.dna_structure.write(file_name,write_json_format=True)

    def write_cando_file(self, file_name):
        """ Write a CanDo .cndo file.
            Arguments:
                file_name (String): The name of the CanDo file to write. 
        """
        cando_writer = CandoWriter(self.dna_structure)
        cando_writer.write(file_name)

    def write_cadnano_file(self, file_name):
        """ Write a caDNAno JSON file.

            Arguments:
                file_name (String): The name of the caDNAno file to write. 
        """
        cadnano_writer = CadnanoWriter(self.dna_structure)
        cadnano_writer.write(file_name)

    def perform_staple_operations(self, staples_arg):
        """ Perform operations on staples.  

            Arguments:
                staples_arg (String): The argument to the staples command-line option.
        """
        tokens = staples_arg.split(",", 1)
        operation = tokens[0]
        retain_staples = []

        # Parse retained staples IDs. 
        if len(tokens) == 2:
            pattern = re.compile('\W')
            retain_tokens = pattern.split(tokens[1])
            if retain_tokens[0] == "retain": 
                retain_colors = [ int(color) for color in retain_tokens[1:] if color != '']
            #__if retain_tokens[0] == "retain"
            retain_staples = self.dna_structure.get_staples_by_color(retain_colors)
        #__if len(tokens) == 2

        # Remove all staple strands except those given in retain_staples[].
        if operation == "delete": 
            self.dna_structure.remove_staples(retain_staples)

        # Generaqte the maximal staple strand set except those given in retain_staples[].
        elif operation == "maximal_set": 
            self.dna_structure.generate_maximal_staple_set(retain_staples)

    #__def perform_staple_operations

    def transform_structure(self, transform):
        """ Apply 3D geometric transformations to a selected set of helices. 

            The format of the transform commands is:
                helices(0,1):rotate(90,0,0),translate(0,0,0);helices(2,3):rotate(0,90,0),translate(0,0,0)
        """
        helices_map = self.dna_structure.structure_helices_map
        self.logger.info("Transform %s" % transform)
        helix_groups = transform.split(";")
        self.logger.info("Number of helix groups %d" % len(helix_groups))

        # Parse helix IDs.
        helix_group_xforms = []
        for helix_group in helix_groups:
            tokens = helix_group.split(":")
            pattern = re.compile(r"[,()]")
            helix_tokens = pattern.split(tokens[0])
            helix_ids = []
            for s in helix_tokens:
                if s == "helices": 
                    continue 
                elif "-" in s:
                    rtoks = s.split("-")
                    start = int(rtoks[0])
                    end = int(rtoks[1])+1
                    for id in xrange(start,end):
                        helix_ids.append(id)
                elif s:
                    helix_ids.append(int(s))
            #__for s in helix_tokens

            # Check helix IDs.
            helices = []
            for hid in helix_ids:
                helix = helices_map.get(hid, None)
                if not helix:
                    self.logger.error("Helix ID %d not found in dna structure." % hid)
                    self.logger.error("DNA Structure has helix IDs %s " % str(helices_map.keys()))
                    return
                helices.append(helix)
                #__if not helix:
            #__for hid in helix_ids
            self.logger.info("Helix group %s" % str(helix_ids))

            # Parse transformations.
            self.logger.info("Transformation \'%s\'" % tokens[1])
            pattern = re.compile(r"[(),]")
            xform_tokens = pattern.split(tokens[1])
            n = 0
            use_connectors = False
            xform = Xform()
            while (n != len(xform_tokens)): 
                s = xform_tokens[n]
                if s == "rotate":
                    rotations = []
                    rotations.append(float(xform_tokens[n+1]))
                    rotations.append(float(xform_tokens[n+2]))
                    rotations.append(float(xform_tokens[n+3]))
                    n += 3
                    xform.add_rotation(rotations)
                    rotations = []
                elif s == "translate":
                    translation = []
                    translation.append(float(xform_tokens[n+1]))
                    translation.append(float(xform_tokens[n+2]))
                    translation.append(float(xform_tokens[n+3]))
                    n += 3
                    xform.set_translation(translation)
                elif s == "connectors":
                    use_connectors = True
                    strand_name = xform_tokens[n+1]
                    n += 1
                #__if s == "rotate"
                n += 1
            #__while (n != len(xform_tokens))

            # Automatically generate the transformation the moves one group of helices to another
            # using the connections of distance crossovers. 
            if use_connectors:
                self.logger.info("Use connectors with strand \'%s\'" % strand_name)
                connector_strands = []
                for strand in self.dna_structure.strands:
                    if strand.is_scaffold:
                        connector_strands.append(strand)
                helix_dist = self.dna_structure.dna_parameters.helix_distance
                xform_from_connectors(connector_strands, helix_ids, helix_dist, xform)
            #__if use_connectors

            helix_group_xforms.append( HelixGroupXform(helices, xform) )
        #__for helix_group in helix_groups

        # Apply the transformation to the dna structure helices.
        apply_helix_xforms(helix_group_xforms) 
    #__def transform_structure

    def set_module_loggers(self, names):
        module_names = names.split(",")
        for module in module_names:
            logger = logging.getLogger(module)
            logger.setLevel(logging.DEBUG)
        #__for module in modules
    #__def set_debugging_loggers

