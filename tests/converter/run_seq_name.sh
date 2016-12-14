#!/bin/bash 
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
 
# This script tests converting a caDNAno design file into an Autodesk Nanodesign
# viewer file. 
#
# The sequence to assign to the design is given by a sequence name.

# Set the cadnano design files directory.
data=../samples

# Set the input caDNAno design file.
fn=fourhelix

# Set the sequence name. 
seq=M13mp18

# Set the sequence name.
# Valid names are: p7308, p7704, p8064, p8100, p8634, M13KO7, p7560 and M13mp18.
seq=M13mp18

# Modify the design with the inserts and deletes given in the caDNAno design file.
# Set to "true" to modify the design.
modify="false"

# Activate debug logging for the given module name(s). 
# Examples: 
#     debug="nanodesign.data.dna_structure,nanodesign.data.dna_structure_helix" - activate debugging for dna_structure and dna_structure_helix modules only
#     debug="nanodesign.data" - activate debugging for all modules under nanodesign.data 
debug=""

outfile=./results/${fn}_viewer.json 

if [ ! -d "results/" ]; then
    mkdir results
fi

../../scripts/converter.py --infile=${data}/${fn}.json  \
                           --informat="cadnano"         \
                           --inseqname=${seq}           \
                           --modify=${modify}           \
                           --debug=${debug}             \
                           --outfile=${outfile}         \
                           --outformat="viewer"

