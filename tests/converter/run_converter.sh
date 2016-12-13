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
# viewer file. No sequence information is created.

# Set the cadnano design file directory.
data=../samples

# Set the input caDNAno design file.
fn=fourhelix

if [ ! -d "results/" ]; then
    mkdir results
fi

# Set the name of the output viewer file.
outfile=./results/${fn}_viewer.json  

../../scripts/converter.py --infile=${data}/${fn}.json \
                           --informat="cadnano"        \
                           --outfile=${outfile}        \
                           --outformat="viewer"

