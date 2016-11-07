# This shell script demonstrates how to execute the visualizer.

# directory path to cadnano design files.
dir=../samples/

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
fn=slottedcross

# set the flag for generating atomic structures.
atomic_model="true"
atomic_model="false"

# execute the visualizer with a cadnano file
fn=slottedcross
fn=fourhelix
fn=42hb
fn=42hb_recode
fn=6hb_recode

if [ $fn == "fourhelix" ]; then
   cmds="helix name=0 rep=geometry color=(1,0.5,0) show=true"
   cmds="strand name=staple_0_26  rep=path  color=(1,0,0) show=true"
   cmds="strand names=start_helices[0]  rep=path  color=(1,0,0) line_width=4.0  show=true"
   cmds="strand names=in_helices[19]  rep=path  color=(1,0,0) line_width=4.0  show=true"
fi

../../scripts/vis.py --infile=${dir}/${fn}.json       \
                     --inseqname=${seqname}           \
                     --atomic_model=${atomic_model}   \
                     --commands="${cmds}" 

