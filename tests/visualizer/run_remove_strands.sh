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
 
#!/bin/sh 
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
 
# This shell script demonstrates how to execute the visualizer to show the deleted staples.
#
# A list of staple colors, matching those defined for staples in the cadnano files, may
# be given to identify the staples retained after deletion.

# Set the directory path to cadnano design files.
dir=../samples/

# Set the sequence name.
seqname=M13mp18

# Set the flag for generating atomic structures.
atomic_model="false"

# Set cadnano file name.
fn=simple

# Set the staples to retain.

if [ $fn == "fourhelix" ]; then
    # The valid staple colors from the fourhelix design file (helix number: [position,color], ...):
    #     1: [8,11184640],[41,243362
    #     0: [26,29184]  ,[60,7536862]
    #     3: [8, 5749504]
    #     2: [46, 243362 ], [60, 16225054] 
   del_cmd="delete,retain=[11184640]"
   del_cmd="delete,retain=[243362]"
   del_cmd="delete,retain=[29184]"
fi

# Set the command to show all stand paths.
cmds="strand name=All  rep=path  show=true"

../../scripts/vis.py --infile=${dir}/${fn}.json       \
                     --inseqname=${seqname}           \
                     --atomic_model=${atomic_model}   \
                     --staples=${del_cmd}             \
                     --commands="${cmds}" 
