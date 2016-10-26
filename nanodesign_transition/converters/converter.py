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
import argparse
from cadnano.reader import CadnanoReader
from cadnano.writer import CadnanoWriter
from cadnano.convert_design import CadnanoConvertDesign
from viewer.writer import ViewerWriter 
from cando.writer import CandoWriter 
from simdna.writer import SimDnaWriter 
from pdbcif.pdb_writer import PdbWriter 
from pdbcif.cif_writer import CifWriter 

try:
    import os.path
    base_path = os.path.abspath( os.path.dirname(__file__) + '/../../' )
    sys.path.append( base_path )
    from nanodesign_transition.dna_structure import DnaStructure
    from nanodesign_transition.parameters import DnaParameters
    from nanodesign_transition.xform import Xform,HelixGroupXform,apply_helix_xforms,xform_from_connectors
    sys.path = sys.path[:-1]
except ImportError as i:
    print "Could not get nanodesign_transition module"
    raise i

from dna_sequence_data import dna_sequence_data

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
        self.logger = None
        self.modify = False
        self.dna_parameters = DnaParameters()

def parse_args():
    """ Parse command-line arguments."""
    parser = argparse.ArgumentParser()
    parser.add_argument("-hd",  "--helixdist",   help="distance between DNA helices")
    parser.add_argument("-if",  "--informat",    help="input file format: cadnano, viewer")
    parser.add_argument("-i",   "--infile",      help="input file")
    parser.add_argument("-is",  "--inseqfile",   help="input sequence file")
    parser.add_argument("-isn", "--inseqname",   help="input sequence name")
    parser.add_argument("-m",   "--modify",      help="create DNA structure using the deleted/inserted bases given in a cadnano design file")
    parser.add_argument("-o",   "--outfile",     help="output file")
    parser.add_argument("-of",  "--outformat",   help="output file format")
    parser.add_argument("-s",   "--staples",     help="staple operations")
    parser.add_argument("-x",   "--transform",   help="apply a transformation to a set of helices")
    return parser.parse_args()

def read_cadnano_file(converter, file_name, seq_file_name, seq_name):
    """ Read in a caDNAno file. 

        Arguments:
            converter (Converter): The converter object that stores various models created when reading from a file.
            file_name (String): The name of the caDNAno file to convert. 
            seq_file_name (String): The name of the CSV file used to assign a DNA base sequence to the DNA structure. 
            seq_name (String): The name of a sequence used to assign a DNA base sequence to the DNA structure.
    """
    cadnano_reader = CadnanoReader()
    converter.cadnano_design = cadnano_reader.read_json(file_name)
    converter.cadnano_convert_design = CadnanoConvertDesign(converter.dna_parameters)
    converter.dna_structure = converter.cadnano_convert_design.create_structure(converter.cadnano_design, converter.modify)

    # Read in staple sequences from a CSV format file.
    if (seq_file_name): 
        _, file_extension = os.path.splitext(seq_file_name)

        if (file_extension == ".csv"): 
            modified_structure = False
            sequence = cadnano_reader.read_csv(seq_file_name)
            converter.cadnano_convert_design.set_sequence(converter.dna_structure, modified_structure, sequence)

    # Assign a sequence using a name.
    if (seq_name): 
        if (seq_name not in dna_sequence_data):
            converter.logger.error("The sequence name %s is not recognized.", seq_name)
        modified_structure = False
        converter.cadnano_convert_design.set_sequence_from_name(converter.dna_structure, modified_structure, seq_name)

def write_viewer_file(converter, file_name):
    """ Write a Nanodesign Viewer file.

        Arguments:
            converter (Converter): The converter object that stores various models created when reading from a file.
            file_name (String): The name of the Nanodesign Viewer file to write. 
    """
    viewer_writer = ViewerWriter(converter.dna_structure, converter.dna_parameters)
    viewer_writer.write(file_name)

def write_pdb_file(converter, file_name):
    """ Write an RCSB PDB-format file.

        Arguments:
            converter (Converter): The converter object that stores various models created when reading from a file.
            file_name (String): The name of the PDB file to write. 
    """
    pdb_writer = PdbWriter(converter.dna_structure)
    pdb_writer.write(file_name)

def write_cif_file(converter, file_name):
    """ Write a RCSB CIF-format file.

        Arguments:
            converter (Converter): The converter object that stores various models created when reading from a file.
            file_name (String): The name of the CIF file to write. 
    """
    cif_writer = CifWriter(converter.dna_structure)
    cif_writer.write(file_name, converter.infile, converter.informat )

def write_simdna_file(converter, file_name):
    """ Write a SimDNA pairs file.

        Arguments:
            converter (Converter): The converter object that stores various models created when reading from a file.
            file_name (String): The name of the SimDNA pairs file to write. 
    """
    simdna_writer = SimDnaWriter(converter.dna_structure)
    simdna_writer.write(file_name)

def write_topology_file(converter, file_name):
    """ Write a DNA topology file.

        Arguments:
            converter (Converter): The converter object that stores various models created when reading from a file.
            file_name (String): The name of the topology file to write. 
    """
    converter.dna_structure.write_topology(file_name, write_json_format=True)

def write_structure_file(converter, file_name):
    """ Write a DNA structure file.
        Arguments:
            converter (Converter): The converter object that stores various models created when reading from a file.
            file_name (String): The name of the structure file to write. 
    """
    converter.dna_structure.write(file_name,write_json_format=True)

def write_cando_file(converter, file_name):
    """ Write a CanDo .cndo file.
        Arguments:
            converter (Converter): The converter object that stores various models created when reading from a file.
            file_name (String): The name of the CanDo file to write. 
    """
    cando_writer = CandoWriter(converter.dna_structure)
    cando_writer.write(file_name)

def write_cadnano_file(converter, file_name):
    """ Write a caDNAno JSON file.

        Arguments:
            converter (Converter): The converter object that stores various models created when reading from a file.
            file_name (String): The name of the caDNAno file to write. 
    """
    cadnano_writer = CadnanoWriter(converter.dna_structure)
    cadnano_writer.write(file_name)

def perform_staple_operations(converter, staples_arg):
    """ Perform operations on staples.  

        Arguments:
            staples_arg (String): The argument to the staples command-line option.
    """
    print("============== perform_staple_operations ============")
    tokens = staples_arg.split(",", 1)
    print(">>> tokens %s" % str(tokens))
    operation = tokens[0]
    retain_staples = []
    
    # Parse retained staples IDs. 
    if len(tokens) == 2:
        pattern = re.compile('\W')
        retain_tokens = pattern.split(tokens[1])
        if retain_tokens[0] == "retain": 
            retain_colors = [ int(color) for color in retain_tokens[1:] if color != '']
        #__if retain_tokens[0] == "retain"
        print(">>> retain_colors %s" % str(retain_colors))
        retain_staples = converter.dna_structure.get_staples_by_color(retain_colors)
        print(">>> Number of retained staples %d" % len(retain_staples))
    #__if len(tokens) == 2

    # Remove all staple strands except those given in retain_staples[].
    if operation == "delete": 
        converter.dna_structure.remove_staples(retain_staples)

#__def perform_staple_operations


def transform_structure(converter, transform):
    """ Apply 3D geometric transformations to a selected set of helices. 

        The format of the transform commands is:
            helices(0,1):rotate(90,0,0),translate(0,0,0);helices(2,3):rotate(0,90,0),translate(0,0,0)
    """
    helices_map = converter.dna_structure.structure_helices_map
    converter.logger.info("Transform %s" % transform)
    helix_groups = transform.split(";")
    converter.logger.info("Number of helix groups %d" % len(helix_groups))

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
                converter.logger.error("Helix ID %d not found in dna structure." % hid)
                converter.logger.error("DNA Structure has helix IDs %s " % str(helices_map.keys()))
                return
            helices.append(helix)
            #__if not helix:
        #__for hid in helix_ids
        converter.logger.info("Helix group %s" % str(helix_ids))

        # Parse transformations.
        converter.logger.info("Transformation \'%s\'" % tokens[1])
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
            converter.logger.info("Use connectors with strand \'%s\'" % strand_name)
            connector_strands = []
            for strand in converter.dna_structure.strands:
                if strand.is_scaffold:
                    connector_strands.append(strand)
            helix_dist = converter.dna_structure.dna_parameters.helix_distance
            xform_from_connectors(connector_strands, helix_ids, helix_dist, xform)
        #__if use_connectors

        helix_group_xforms.append( HelixGroupXform(helices, xform) )
    #__for helix_group in helix_groups

    # Apply the transformation to the dna structure helices.
    apply_helix_xforms(helix_group_xforms) 
#__def transform_structure

def _setup_logging():
    """ Set up logging."""
    logger = logging.getLogger('nanodesign.converter')
    logger.setLevel(logging.INFO)

    # create console handler and set format
    console_handler = logging.StreamHandler()
    formatter = logging.Formatter('[%(name)s] %(levelname)s - %(message)s')
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)
    return logger

converter_read_map = { ConverterFileFormats.CADNANO    : read_cadnano_file }

converter_write_map = { ConverterFileFormats.VIEWER    : write_viewer_file,
                        ConverterFileFormats.CADNANO   : write_cadnano_file,
                        ConverterFileFormats.CANDO     : write_cando_file,
                        ConverterFileFormats.CIF       : write_cif_file,
                        ConverterFileFormats.PDB       : write_pdb_file,
                        ConverterFileFormats.SIMDNA    : write_simdna_file,
                        ConverterFileFormats.STRUCTURE : write_structure_file,
                        ConverterFileFormats.TOPOLOGY  : write_topology_file
                      }

def main():
    logger = _setup_logging()
    converter = Converter()
    converter.logger = logger
    infiles = []

    # Process command-line arguments.
    args = parse_args()

    if args.infile == None:
        logger.error("No input file name given.")
    else:
        logger.info("Input file name: %s" % args.infile)
        converter.infile = args.infile

    if args.informat == None:
        logger.error("No input file format given.")
    elif (args.informat not in ConverterFileFormats.names):
        logger.error("Unknown input file format given: %s" % args.informat)
    else:
        logger.info("Input file format: %s" % args.informat)
        converter.informat = args.informat

    if args.modify:
        logger.info("Create a DNA structure using deleted/inserted bases from the caDNAno design file.")
        converter.modify = (args.modify.lower() == "true")

    if args.helixdist:
        converter.dna_parameters.helix_distance = float(args.helixdist)
        logger.info("Set the distance between adjacent helices to %g" % converter.dna_parameters.helix_distance)

    if args.outfile == None:
        logger.error("No output file name given.")
    else:
        logger.info("Output file name: %s" % args.infile)

    if args.outformat == None:
        logger.error("No output file format given.")
    elif (args.outformat not in  ConverterFileFormats.names):
        logger.error("Unknown output file format given: \'%s\'" % args.outformat)
    else:
        logger.info("Output file format: %s" % args.outformat)
        # Make the helix distance a bit larger to better visualization.
        if args.outformat == ConverterFileFormats.VIEWER:
            converter.dna_parameters.helix_distance = 2.50

    # read the input file
    converter_read_map[args.informat](converter, args.infile, args.inseqfile, args.inseqname)

    # perform staple operations (e.g., delete, generate maximal set, etc.) 
    if args.staples:
        perform_staple_operations(converter, args.staples)

    # appy a 3D transformation to the geometry of selected helices.
    if args.transform:
        transform_structure(converter, args.transform)

    # write the output file
    converter_write_map[args.outformat](converter, args.outfile)

if __name__ == '__main__':
    main()

