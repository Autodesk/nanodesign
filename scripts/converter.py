#!/usr/bin/env python

import logging
import argparse
import os.path

try:
    from nanodesign.converters.converter import Converter,ConverterFileFormats
except ImportError:
    import sys
    base_path = os.path.abspath( os.path.join( os.path.dirname(os.path.abspath( __file__)), '../'))
    sys.path.append(base_path)
    from nanodesign.converters.converter import Converter,ConverterFileFormats
    sys.path = sys.path[:-1]

# Define the map between file formats and the functions that read files in that format. 
converter_read_map = { ConverterFileFormats.CADNANO    : 'read_cadnano_file' }

# Define the map between file formats and the functions that write files in that format. 
converter_write_map = { ConverterFileFormats.VIEWER    : 'write_viewer_file',
                        ConverterFileFormats.CADNANO   : 'write_cadnano_file',
                        ConverterFileFormats.CANDO     : 'write_cando_file',
                        ConverterFileFormats.CIF       : 'write_cif_file',
                        ConverterFileFormats.PDB       : 'write_pdb_file',
                        ConverterFileFormats.SIMDNA    : 'write_simdna_file',
                        ConverterFileFormats.STRUCTURE : 'write_structure_file',
                        ConverterFileFormats.TOPOLOGY  : 'write_topology_file'
                      }

def parse_args():
    """ Parse command-line arguments."""
    parser = argparse.ArgumentParser()
    parser.add_argument("-dbg", "--debug",       help="set modules debugging logger")
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

def main():
    logger = logging.getLogger('nanodesign.converter')
    converter = Converter()
    infiles = []

    # Process command-line arguments.
    args = parse_args()

    if args.debug:
        converter.set_module_loggers(args.debug) 

    if args.infile == None:
        logger.error("No input file name given.")
    else:
        logger.info("Input file name %s" % args.infile)
        converter.infile = args.infile

    if args.informat == None:
        logger.error("No input file format given.")
    elif (args.informat not in ConverterFileFormats.names):
        logger.error("Unknown input file format given: %s" % args.informat)
    else:
        logger.info("Input file format %s" % args.informat)
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
        logger.info("Output file name %s" % args.outfile)

    if args.outformat == None:
        logger.error("No output file format given.")
    elif (args.outformat not in  ConverterFileFormats.names):
        logger.error("Unknown output file format given \'%s\'" % args.outformat)
    else:
        logger.info("Output file format %s" % args.outformat)
        # Make the helix distance a bit larger to better visualization.
        if args.outformat == ConverterFileFormats.VIEWER:
            converter.dna_parameters.helix_distance = 2.50

    # read the input file
    read_function = getattr( converter, converter_read_map[args.informat] )
    read_function(args.infile, args.inseqfile, args.inseqname)

    # perform staple operations (e.g., delete, generate maximal set, etc.) 
    if args.staples:
        converter.perform_staple_operations(args.staples)

    # appy a 3D transformation to the geometry of selected helices.
    if args.transform:
        converter.transform_structure(args.transform)

    # write the output file
    write_function = getattr( converter, converter_write_map[args.outformat] )
    write_function(args.outfile)

if __name__ == '__main__':
    main()
