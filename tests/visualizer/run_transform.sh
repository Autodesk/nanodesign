# This shell script demonstrates how to execute the visualizer and show a transformed structure.
# The geometry of the helices in the design are rotated and translated.

# Set the directory path to cadnano design files.
dir=../samples/

# Set the sequence name.
seqname=M13mp18

# Set cadnano file name
fn=fourhelix

# Set the flag for generating atomic structures.
atomic_model="false"

# Set the helix IDs and transformations.
#
# The 'helices' option sepcifies a list of caDNAno helix IDs .
#
# The transformation is given by
#
#    rotate(rx,ry,rz) - rotates about x, y and z axes.
#
#    translate(tx,ty,tz) - translates by tx, ty and tz.

if [ $fn == "fourhelix" ]; then 
    transform="helices(0,1):rotate(90,0,0),translate(0.5,0,0);helices(2,3):rotate(0,90,0),translate(0,0,0)"
    cmds="helix name=All  rep=geometry  show=true"
fi

../../scripts/vis.py --infile=${dir}/${fn}.json       \
                     --inseqname=${seqname}           \
                     --atomic_model=${atomic_model}   \
                     --transform=${transform}         \
                     --commands="${cmds}" 

