# This shell script demonstrates how to execute the visualizer.

# directory path to cadnano design files.
dir=../../tests/samples/

# sequence file name
seqfile=fourhelix_M13mp18

# sequence name
seqname=M13mp18

# various cadnano file names
fn=42hb
fn=AutodeskA
fn=doublegear
fn=fourhelix
fn=hc-test-1
fn=nanorobot.v2
fn=Nature09_squarenut_no_joins
fn=robot_v1.9_bent_2
fn=Science09_beachball_v1
fn=slottedcross-1

# set the flag for generating atomic structures.
atomic_model="true"

# execute the visualizer with a cadnano file
fn=fourhelix
vis.py --infile=${dir}/${fn}.json  --inseqname=${seqname}  --atomic_model=${atomic_model}

