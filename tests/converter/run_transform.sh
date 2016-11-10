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

