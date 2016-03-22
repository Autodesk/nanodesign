#!/usr/bin/env python
""" 
This module is used to convert DNA design files into other file formats.
"""
import os
import sys
import json
import logging
import argparse
from cadnano.reader import CadnanoReader
from cadnano.convert_design import CadnanoConvertDesign
from viewer.writer import ViewerWriter 
from cando.writer import CandoWriter 
from nanodesign.dna_structure import DnaStructure
from dna_sequence_data import dna_sequence_data

class ConverterFileFormats(object):
    """ File format names to convert to/from. """
    UNKNOWN = "unknown"
    CADNANO = "cadnano"
    CANDO   = "cando"
    VIEWER  = "viewer"
    names = [ CADNANO, CANDO, VIEWER ]

class Converter(object):
    """ This class stores objects for various models created when reading from a file.
    """
    def __init__(self):
        self.cadnano_design = None 
        self.dna_structure = None 
        self.cadnano_convert_design = None
        self.logger = None

def parse_args():
    """ Parse command-line arguments."""
    parser = argparse.ArgumentParser()
    parser.add_argument("-if", "--informat",   help="input file format: cadnano, viewer")
    parser.add_argument("-i",  "--infile",     help="input file")
    parser.add_argument("-is", "--inseqfile",  help="input sequence file")
    parser.add_argument("-isn","--inseqname",  help="input sequence name")
    parser.add_argument("-o",  "--outfile",    help="output file")
    parser.add_argument("-of", "--outformat",  help="output file format")
    return parser.parse_args()

def read_cadnano_file(converter, file_name, seq_file_name, seq_name):
    """ Read in a caDNAno file."""
    cadnano_reader = CadnanoReader()
    converter.cadnano_design = cadnano_reader.read_json(file_name)
    converter.cadnano_convert_design = CadnanoConvertDesign()
    converter.dna_structure = converter.cadnano_convert_design.create_structure(converter.cadnano_design)

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
    viewer_writer = ViewerWriter(converter.dna_structure)
    viewer_writer.write(file_name)

def write_cando_file(converter, file_name):
    """ Write a CanDo .cndo file."""
    cando_writer = CandoWriter(converter.dna_structure)
    cando_writer.write(file_name)

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

converter_read_map = { ConverterFileFormats.CADNANO : read_cadnano_file }
converter_write_map = { ConverterFileFormats.VIEWER : write_viewer_file,
                        ConverterFileFormats.CANDO  : write_cando_file
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

    if args.informat == None:
        logger.error("No input file format given.")
    elif (args.informat not in ConverterFileFormats.names):
        logger.error("Unknown input file format given: %s" % args.informat)
    else:
        logger.info("Input file format: %s" % args.informat)

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

    # read the input file
    converter_read_map[args.informat](converter, args.infile, args.inseqfile, args.inseqname)

    # write the output file
    converter_write_map[args.outformat](converter, args.outfile)

if __name__ == '__main__':
    main()

