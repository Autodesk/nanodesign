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
 
# This shell script demonstrates how to execute the visualizer with command-line commands.

# Set the directory path to cadnano design files.
dir=../samples/

# Set the sequence name.
seqname=M13mp18

# Set the cadnano file name.
fn=fourhelix

# Set the flag for generating atomic structures.
atomic_model="true"

# Examples of some visualizer commands.
if [ $fn == "fourhelix" ]; then
   cmds="helix name=0 rep=geometry color=(1,0.5,0) show=true"
   cmds="strand name=staple_0_26  rep=path  color=(1,0,0) show=true"
   cmds="strand names=start_helices[0]  rep=path  color=(1,0,0) line_width=4.0  show=true"
   cmds="strand names=in_helices[19]  rep=path  color=(1,0,0) line_width=4.0  show=true"
fi

../../scripts/vis.py --infile=${dir}/${fn}.json     \
                     --inseqname=${seqname}         \
                     --atomic_model=${atomic_model} \
                     --commands="${cmds}"

