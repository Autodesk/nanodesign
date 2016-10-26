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
fn=slottedcross-1
#dir=~/
#fn=C_G_1W_fixed

# set the flag for generating atomic structures.
atomic_model="true"
atomic_model="false"

# execute the visualizer with a cadnano file
fn=protractor_30_98_v4
fn=fourhelix_inserts_1
fn=fourhelix_inserts_2
fn=fourhelix

#cmds="model rep=Geometry show=true"
#cmds="helix name=2  rep=base_positions  show=true;helix name=3  rep=geometry  show=true;helix name=2  rep=geometry  show=true"
#cmds="helix name=1  rep=base_positions  show=true;helix name=1  rep=geometry  show=true;helix name=0  rep=geometry  show=true"
cmds="strand name=Scaffold_55_7  rep=path  show=true;strand name=Scaffold_55_7  rep=connectors  show=true "

transform="helices(0,1):rotate(90,0,0),translate(0,0,0)"
transform="helices(0,1):rotate(0,90,0),translate(0,0,0)"
transform="helices(0,1):rotate(0,0,90),translate(0,0,0)"
transform="helices(2,3):rotate(0,0,90),translate(0,0,0)"
transform="helices(0,1):rotate(0,0,90),translate(0,0,0)"
transform=""

fn=slottedcross
transform="helices(0-41):rotate(180,0,90),translate(0,0,0)"
transform="helices(0-41):rotate(0,90,0),translate(0,0,0)"
transform="helices(0-41):rotate(0,180,90),translate(0,0,40)"
transform="helices(0-41):rotate(0,0,0),translate(0,0,30)"
#transform=""
#transform="helices(0-41):connectors(scaffold)"
transform="helices(0-41):rotate(90,0,90),translate(0,0,41.4)"
transform=""

../../scripts/vis.py --infile=${dir}/${fn}.json       \
                     --inseqname=${seqname}           \
                     --atomic_model=${atomic_model}   \
                     --transform=${transform}         \
                     --commands="${cmds}" 

