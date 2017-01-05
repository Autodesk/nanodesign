#!/usr/bin/env python

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

"""This module demonstrates using the Nanodesign interface to calculate some statistics
   for the strands in a design.

    Using the length of each strand the number of instances of each strand length 
    are counted. The counts for each length are then sorted and printed.

    The number of helices each strand visits is then calculated using its list of bases.
    The counts for each helix are then sorted and printed.
"""
from collections import OrderedDict
import logging
import os
import re
import sys

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

from nanodesign.converters import Converter

def read_file(file_name, seq_name):
    """ Read in a cadnano file. """
    converter = Converter()
    seq_file = None
    converter.read_cadnano_file(file_name, seq_file, seq_name)
    return converter

def main():
    # Set caDNAno file name.
    file_name = "../tests/samples/Nature09_squarenut_no_joins.json"

    # Set sequence to assign to scaffold.
    seq_name = "M13mp18"

    # Read cadnano file and create dna structure.
    converter = read_file(file_name, seq_name)
    dna_structure = converter.dna_structure

    # Count the number of instances of each strand length.
    strand_lengths = OrderedDict()
    for strand in dna_structure.strands:
        num_bases = len(strand.tour)
        if num_bases not in strand_lengths:
            strand_lengths[num_bases] = 0
        strand_lengths[num_bases] += 1
    #__for strand in self.strands
    print("\nStrand length counts:")
    for length, count in sorted(strand_lengths.iteritems(), key=lambda (k,v): (v,k)):
        print 'Length {:>4}  Count {:>4}'.format(length, count)
    #__for length in strand_lengths

    # Calculate the number of helices each strand visits.
    strand_helix = OrderedDict()
    for strand in dna_structure.strands:
        id = strand.id
        helix_ids = set()
        for base in strand.tour:
            helix_ids.add(base.h)
        num_helices = len(helix_ids)
        if num_helices not in strand_helix:
            strand_helix[num_helices] = 0
        strand_helix[num_helices] += 1
    #__for strand in dna_structure.strands
    print("\nStrand helix counts:")
    for num_helices, count in sorted(strand_helix.iteritems(), key=lambda (k,v): (v,k)):
        print 'Number of helices {:>4}  Count {:>4}'.format(num_helices, count) 
    #__for num_helices in strand_helix

if __name__ == '__main__':
    main()

