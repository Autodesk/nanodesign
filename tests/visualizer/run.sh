# This shell script demonstrates how to execute the visualizer.

# Set the directory path to cadnano design files.
dir=../samples/

# Set the sequence name.
seqname=M13mp18

# Set the input cadnano file name.
fn=fourhelix

# Set the flag for generating atomic structures.
atomic_model="true"

../../scripts/vis.py --infile=${dir}/${fn}.json       \
                     --inseqname=${seqname}           \
                     --atomic_model=${atomic_model} 

