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
from pdb.pdb_writer import PdbWriter 
from pdb.cif_writer import CifWriter 

try:
    import os.path
    base_path = os.path.abspath( os.path.dirname(__file__) + '/../../' )
    sys.path.append( base_path )
    from nanodesign_transition.dna_structure import DnaStructure
    from nanodesign_transition.parameters import DnaParameters
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
    names = [ CADNANO, CANDO, CIF, PDB, STRUCTURE, TOPOLOGY, VIEWER ]

class Converter(object):
    """ This class stores objects for various models created when reading from a file.
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
        self.helix_distance = DnaParameters.helix_distance

def parse_args():
    """ Parse command-line arguments."""
    parser = argparse.ArgumentParser()
    parser.add_argument("-hd",  "--helixdist",  help="distance between DNA helices")
    parser.add_argument("-if",  "--informat",   help="input file format: cadnano, viewer")
    parser.add_argument("-i",   "--infile",     help="input file")
    parser.add_argument("-is",  "--inseqfile",  help="input sequence file")
    parser.add_argument("-isn", "--inseqname",  help="input sequence name")
    parser.add_argument("-m",   "--modify",     help="create DNA structure using the deleted/inserted bases given in a cadnano design file")
    parser.add_argument("-o",   "--outfile",    help="output file")
    parser.add_argument("-of",  "--outformat",  help="output file format")
    return parser.parse_args()

def read_cadnano_file(converter, file_name, seq_file_name, seq_name):
    """ Read in a caDNAno file. """
    cadnano_reader = CadnanoReader()
    converter.cadnano_design = cadnano_reader.read_json(file_name)
    converter.cadnano_convert_design = CadnanoConvertDesign()
    converter.dna_structure = converter.cadnano_convert_design.create_structure(converter.cadnano_design, converter.modify,
        converter.helix_distance)

    if (seq_file_name): 
        _, file_extension = os.path.splitext(seq_file_name)

        if (file_extension == ".csv"): 
            modified_structure = False
            sequence = cadnano_reader.read_csv(seq_file_name)
            converter.cadnano_convert_design.set_sequence(modified_structure, sequence)

    if (seq_name): 
        if (seq_name not in dna_sequence_data):
            converter.logger.error("The sequence name %s is not recognized.", seq_name)
        modified_structure = False
        converter.cadnano_convert_design.set_sequence_from_name(modified_structure, seq_name)

def write_viewer_file(converter, file_name):
    """ Write a DNA Design viewer file."""
    viewer_writer = ViewerWriter(converter.dna_structure, converter.helix_distance)
    viewer_writer.write(file_name)

def write_pdb_file(converter, file_name):
    """ Write a PDB file."""
    pdb_writer = PdbWriter(converter.dna_structure)
    pdb_writer.write(file_name)

def write_cif_file(converter, file_name):
    """ Write a CIF file."""
    cif_writer = CifWriter(converter.dna_structure)
    cif_writer.write(file_name, converter.infile, converter.informat )

def write_simdna_file(converter, file_name):
    """ Write a SimDNA pairs file."""
    simdna_writer = SimDnaWriter(converter.dna_structure)
    simdna_writer.write(file_name)

def write_topology_file(converter, file_name):
    """ Write a DNA topology file."""
    converter.dna_structure.write_topology(file_name,write_json_format=True)

def write_structure_file(converter, file_name):
    """ Write a DNA structure file."""
    converter.dna_structure.write(file_name,write_json_format=True)

def write_cando_file(converter, file_name):
    """ Write a CanDo .cndo file."""
    cando_writer = CandoWriter(converter.dna_structure)
    cando_writer.write(file_name)

def write_cadnano_file(converter, file_name):
    """ Write a caDNAno JSON file."""
    cadnano_writer = CadnanoWriter(converter.dna_structure)
    cadnano_writer.write(file_name)

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
        converter.helix_distance = float(args.helixdist)
        logger.info("Set the distance between adjacent helices to %g" % converter.helix_distance)

    if args.outfile == None:
        logger.error("No output file name given.")
    else:
        logger.info("Output file name: %s" % args.infile)

    if args.outformat == None:
        logger.error("No output file format given.")
    elif (args.outformat not in  ConverterFileFormats.names):
        logger.error("Unknown output file format given: %s" % args.outformat)
    else:
        logger.info("Output file format: %s" % args.outformat)
        # Make the helix distance a bit larger to better visualization.
        if args.outformat == ConverterFileFormats.VIEWER:
            converter.helix_distance = 2.50
    # read the input file
    converter_read_map[args.informat](converter, args.infile, args.inseqfile, args.inseqname)

    # write the output file
    converter_write_map[args.outformat](converter, args.outfile)

if __name__ == '__main__':
    main()

