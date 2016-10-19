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
fn=convex_triangle
fn=doublegear
fn=fourhelix
fn=fourhelix_deletes
fn=fourhelix_deletes_1
fn=fourhelix_inserts
fn=fourhelix_inserts_1
fn=fourhelix_ss_inserts
fn=hc-test-1
fn=nanorobot.v2
fn=Nature09_squarenut_no_joins
fn=protractor_30_98_v4
fn=railedbridge
fn=robot_v1.9_bent_2
fn=Science09_beachball_v1
fn=slottedcross-1
#dir=~/
#fn=C_G_1W_fixed

# set the flag for generating atomic structures.
atomic_model="false"
atomic_model="true"

# execute the visualizer with a cadnano file
fn=fourhelix
fn=protractor_30_98_v4
fn=fourhelix_inserts_1
fn=fourhelix_inserts_2
#cmds="model rep=Geometry show=true"
vis.py --infile=${dir}/${fn}.json  --inseqname=${seqname}  --atomic_model=${atomic_model} --commands="${cmds}" 

