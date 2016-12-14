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
# The geometry of the helices in the design are rotated and translated.

# Set the cadnano design files directory.
data=../samples

# Set the input caDNAno design file.
fn=fourhelix

# Set the helix IDs and transformations. 
#
# The 'helices' option sepcifies a list of caDNAno helix IDs .
#
# The transformation is given by 
#
#    rotate(rx,ry,rz) - rotates about x, y and z axes.
#
#    translate(tx,ty,tz) - translates by tx, ty and tz.
transform="helices(0,1):rotate(90,0,0),translate(0,0,0)"
transform="helices(0,1):rotate(90,0,0),translate(0.5,0,0);helices(2,3):rotate(0,90,0),translate(0,0,0)"
transform="helices(1):rotate(0,0,90),translate(0,0,0)"
transform="helices(0,1):rotate(0,0,90),translate(0,0,0)"
transform="helices(0,1-3):rotate(0,0,90),translate(0,0,0)"

outfile=./results/${fn}_viewer.json  

if [ ! -d "results/" ]; then
    mkdir results
fi

../../scripts/converter.py --infile=${data}/${fn}.json \
             --informat="cadnano"        \
             --transform=${transform}    \
             --outfile=${outfile}        \
             --outformat="viewer"

