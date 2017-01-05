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

"""This module demonstrates using the Nanodesign interface to search the domains created from a caDNAno
   design file for a given sequence.

   Each domain in the design is printed with its seqeunce. If the domain contains the search sequence then
   the location of the first occurence of the search sequence is highlighed by a line above it.
"""
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
    file_name = "../tests/samples/fourhelix.json"

    # Set sequence to assign to scaffold.
    seq_name = "M13mp18"

    # Set the search sequence.
    search_seq = "GGT"
  
    # Read cadnano file and create dna structure.
    converter = read_file(file_name, seq_name)
    dna_structure = converter.dna_structure

    # Compute auxillary data to calculate domains.
    dna_structure.compute_aux_data()

    # Search domains.
    print("\nSearch for sequence %s: " % search_seq)
    pattern = re.compile(search_seq)
    for domain in dna_structure.domain_list:
        id = domain.id
        seq = domain.sequence
        match = pattern.search(seq)
        dstr = "Domain ID %4d  Seqence %s" % (id, seq)
        if match:
            loc = [ "_" if i >= match.start() and i < match.end() else " " for i in xrange(0,len(seq)) ] 
            lstr = "".join(loc)
            print '{:>{width}}'.format(lstr, width=len(dstr))
        print '{0}'.format(dstr)
        #__if match
    #__for domain in dna_structure.domain_list

if __name__ == '__main__':
    main()

