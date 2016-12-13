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
 
# This script tests running the converter with commands to perform staple operations: 
#
#     - delete staples 
#
#     - generate a maximal staple set 
#
# for a caDNAno design. A caDNAno design file is written with the modified staples.
#
# A list of staple colors, matching those defined for staples in the cadnano files, may 
# be given to identify the staples retained after an operation.

# Set the cadnano design directory.
data=../samples

# Set the input caDNAno design file.
fn=fourhelix

# Set the staple operation.
op=delete
op=maximal_set

# Set staple commands depending on the input file name.
if [ $fn == "fourhelix" ]; then
    # Delete staples.
    if [ $op == "delete" ]; then
        del_cmd="delete,retain=[243362]"
        del_cmd="delete"
        del_cmd="delete,retain=[243362,29184,7536862,5749504,16225054]"
        del_cmd="delete,retain=[11184640]"
        del_cmd="delete,retain=[29184]"
        del_cmd="delete,retain=[11184640,243362]"
        staple_cmd=${del_cmd}
        outfile=./results/${fn}_nd_del.json  
    fi

    # Create maximal staple set.
    if [ $op == "maximal_set" ]; then
        mset_cmd="maximal_set,retain=[11184640]"
        mset_cmd="maximal_set,retain=[243362]"
        mset_cmd="maximal_set,retain=[29184]"
        mset_cmd="maximal_set,retain=[11184640,243362]"
        mset_cmd="maximal_set"
        staple_cmd=${mset_cmd}
        outfile=./results/${fn}_nd_max.json  
    fi
else 
    if [ $op == "maximal_set" ]; then
        mset_cmd="maximal_set"
        staple_cmd=${mset_cmd}
        outfile=./results/${fn}_nd_max.json  
    fi

    if [ $op == "delete" ]; then
        del_cmd="delete"
        staple_cmd=${del_cmd}
        outfile=./results/${fn}_nd_del.json
    fi
fi

if [ ! -d "results/" ]; then
    mkdir results
fi

../../scripts/converter.py --infile=${data}/${fn}.json   \
                           --informat="cadnano"          \
                           --staples=${staple_cmd}       \
                           --outfile=${outfile}          \
                           --outformat="cadnano"

