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
import re
import sys
import json
import logging
import argparse


try:
    import nanodesign
except ImportError:
    import sys
    base_path = os.path.abspath( os.path.join( os.path.dirname(os.path.abspath( __file__)), '../'))
    sys.path.append(base_path)
    import nanodesign
    # if it fails now, we let the exception go all the way up to halt execution.
    # TODO (JMS 10/26/16): add better reporting of the import error.
    sys.path = sys.path[:-1]

from nanodesign.visualizer.model import VisModel
from nanodesign.converters import Converter
from nanodesign.converters.pdbcif.atomic_structure import AtomicStructure  

def parse_args():
    """ Parse command-line arguments. """
    parser = argparse.ArgumentParser()
    parser.add_argument("-a",  "--atomic_model",   help="generate atomic model")
    parser.add_argument("-c",  "--commands",       help="commands")
    parser.add_argument("-cf", "--cmdfile",        help="command file")
    parser.add_argument("-i",  "--infile",         help="input file")
    parser.add_argument("-is", "--inseqfile",      help="input sequence file")
    parser.add_argument("-isn","--inseqname",      help="input sequence name")
    parser.add_argument("-x",   "--transform",  help="apply a transformation to a set of helices")
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

def transform_structure(logger, dna_structure, transform):
    """ Apply 3D geometric transformations to a selected set of helices. 

        The format of the transform commands is:
            helices(0,1):rotate(90,0,0),translate(0,0,0);helices(2,3):rotate(0,90,0),translate(0,0,0)

        TODO: This code is duplicated in converter.py. Change it to use that code.
    """
    logger.info("Transform %s" % transform)
    helix_groups = transform.split(";")
    logger.info("Number of helix groups %d" % len(helix_groups))
    helix_group_xforms = []
    helices_map = dna_structure.structure_helices_map
    # Parse helix IDs.
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
        logger.info("Helix group %s" % str(helix_ids))

        # Check helix IDs.
        helices = []
        for hid in helix_ids:
            helix = helices_map.get(hid, None)
            if not helix:
                logger.error("Helix ID %d not found in dna structure." % hid)
                logger.error("DNA Structure has helix IDs %s " % str(helices_map.keys()))
                return
            helices.append(helix)
            #__if not helix:
        #__for hid in helix_ids

        # Parse transformations.
        logger.info("Transformation \'%s\'" % tokens[1])
        pattern = re.compile(r"[(),]")
        xform_tokens = pattern.split(tokens[1])
        xform = Xform()
        n = 0
        use_connectors = False
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

        # Create a transformation from connectors.
        if use_connectors:
            logger.info("Use connectors with strand \'%s\'" % strand_name)
            connector_strands = []
            for strand in dna_structure.strands:
                if strand.is_scaffold:
                    connector_strands.append(strand)
            helix_dist = dna_structure.dna_parameters.helix_distance
            xform_from_connectors(connector_strands, helix_ids, helix_dist, xform)
        #__if use_connectors

        helix_group_xforms.append( HelixGroupXform(helices, xform) )
    #__for helix_group in helix_groups

    # Apply the transformation to the dna structure helices.
    apply_helix_xforms(helix_group_xforms) 
#__def transform_structure

def read_file(args, logger):
    """ Read in a cadnano file. """
    converter = Converter()
    converter.logger = logger

    if args.infile == None:
        logger.error("No input file name given.")
        sys.exit(1)
    else:
        logger.info("Input file name: %s" % args.infile)
        infile = args.infile


    converter.read_cadnano_file( args.infile, args.inseqfile, args.inseqname )

    return converter

def main():
    logger = _setup_logging()

    # Get command-line arguments.
    args = parse_args()

    # Read cadnano file and create dna structure.
    converter = read_file(args, logger)

    # appy a 3D transformation to the geometry of selected helices.
    if args.transform:
        converter.transform_structure(args.transform)

    dna_structure = converter.dna_structure
    # No longer need the converter.

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

