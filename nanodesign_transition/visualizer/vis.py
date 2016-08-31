#!/usr/bin/env python
""" 
This module is used to visualize a 3D DNA structure defined by a caDNAno design file.

Virtual helix and strand information is read from a caDNAno design file and used to generate 3D coordinates of virtual 
helices and their base positions based on the design lattice type. These coordinates are then used to visualize the
geometry of virtual helices as cylinders, scaffold and staple strands as they wind through the virtual helices, base 
coordinates and coordinate frames, and domains calculated from base pairing. The backbone P atoms and atom bonds of the 
atomic model generated for the DNA structure can also be visualized.

The 3D geometry of the DNA structure can be visualized by virtual helix or strand basis by selecting an ID from a popup menu.
Virtual helices are selected from a menu using their caDNAno numbering; strands are selected using a name derived from their
type (scaffold or staple), the virtual helix number they start from and the base position within the virtual helix they start 
from. Strand names are sorted by virtual helix number. 

Selecting an ID from a popup menu displays its geometry in the grapics window. Items selected from a menu display a '+' after 
it to indicate it has been selected and is visible. Selecting an item with a '+' will then deselect it and remove it from the 
graphics window. All menus have a 'All' and 'None' entries that are used to select all or none of the items in the menu, 
respectively.

Geometry displayed in the graphics window can be interactively rotated, translated and scaled using the mouse. Geometry can 
also be interactively selected using the mouse. The selection is a 3D operation that highlights and prints any meaningful 
component of the geometry. For example, selecting a bond on a atom bonds geometry of the atomic model highlights all the 
bonds belonging to the base and prints its strand name, sequence number and base name.

Each menu operation generates a visualization command that is written to a file named 'vis.cmd'. These commands can be saved 
to a file and read in to perform operations when the visualizer starts. Commands can also be executed from the command line 
as a semicolon-separated string.
"""
import os
import sys
import json
import logging
import argparse
from model import VisModel

# Set path and import reader and converter.
base_path = os.path.abspath( os.path.dirname(__file__) + '/../../nanodesign_transition/converters' )
sys.path.append( base_path )
from cadnano.reader import CadnanoReader
from cadnano.convert_design import CadnanoConvertDesign
from dna_sequence_data import dna_sequence_data
sys.path = sys.path[:-1]

# Set path and import atomic structure.
base_path = os.path.abspath( os.path.dirname(__file__) + '/../../' )
sys.path.append( base_path )
from nanodesign_transition.atomic_structure import AtomicStructure  
sys.path = sys.path[:-1]

def parse_args():
    """ Parse command-line arguments. """
    parser = argparse.ArgumentParser()
    parser.add_argument("-a",  "--atomic_model",   help="generate atomic model")
    parser.add_argument("-c",  "--commands",       help="commands")
    parser.add_argument("-cf", "--cmdfile",        help="command file")
    parser.add_argument("-i",  "--infile",         help="input file")
    parser.add_argument("-is", "--inseqfile",      help="input sequence file")
    parser.add_argument("-isn","--inseqname",      help="input sequence name")
    return parser.parse_args()

def _setup_logging():
    """ Set up logging. """
    logger = logging.getLogger("visualizer")
    logger.setLevel(logging.INFO)

    # Create console handler and set format.
    console_handler = logging.StreamHandler()
    formatter = logging.Formatter('[%(name)s] %(levelname)s - %(message)s')
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)
    return logger

def read_file(args, logger):
    """ Read in a cadnano file. """

    if args.infile == None:
        logger.error("No input file name given.")
        sys.exit(1)
    else:
        logger.info("Input file name: %s" % args.infile)
        infile = args.infile

    # Read the cadnano file.
    cadnano_reader = CadnanoReader()
    cadnano_design = cadnano_reader.read_json(args.infile)

    # Convert the design.
    convert_design = CadnanoConvertDesign()
    dna_structure = convert_design.create_structure(cadnano_design)

    # Read a sequence file.
    if args.inseqfile:
        logger.info("Input sequence file name: %s" % args.inseqfile)
        _, file_extension = os.path.splitext(args.inseqfile)

        if (file_extension == ".csv"):
            modified_structure = False
            sequence = cadnano_reader.read_csv(args.inseqfile)
            convert_design.set_sequence(modified_structure, sequence)
    #__if (seq_file_name)

    elif (args.inseqname):
        if (args.inseqname not in dna_sequence_data):
            logger.error("The sequence name %s is not recognized.", args.inseqname)
        modified_structure = False
        convert_design.set_sequence_from_name(modified_structure, args.inseqname)
    return dna_structure 

def main():
    logger = _setup_logging()

    # Get command-line arguments.
    args = parse_args()

    # Read cadnano file and create dna structure.
    dna_structure = read_file(args, logger)

    # Generate atomic models of the dna structure.
    if not args.atomic_model or (args.atomic_model == "true"):
        atomic_structure = AtomicStructure(dna_structure)
    else:
        atomic_structure = None

    # Initialize visualization.
    vis_model = VisModel(args.infile, args.cmdfile, args.commands, dna_structure, atomic_structure)

    # Start interactive visualization. 
    vis_model.start_interactive()

if __name__ == '__main__':
    main()

