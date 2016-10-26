#-------------------------------------------------------------------------------------#
#                     run the converter with a sequence name.                         #
#-------------------------------------------------------------------------------------#

# Set the cadnano design files directory.
data=../samples

# A list of various cadnano design file names.
fn=42hb
fn=aNANO_3D_7_14_final
fn=fourhelix
fn=fourhelix_deletes
fn=icosahedron
fn=protractor_30_98_v4
fn=robot_v1.9_bent_2
fn=slottedcross-1
fn=slottedcross

# A list of various sequence names.
seq=p7308
seq=p7704
seq=p8064
seq=p8100
seq=p8634
seq=M13KO7
seq=p7560
seq=M13mp18

# Set to true to modify structure with inserts and deletes.
#modify="false"
modify="true"

# Execute the converter.
fn=fourhelix
fn=slottedcross

outfile=./results/${fn}_viewer.json 

if [ ! -d "results/" ]; then
    mkdir results
fi

../../scripts/converter.py --infile=${data}/${fn}.json      \
             --informat="cadnano"             \
             --inseqname=${seq}               \
             --modify=${modify}               \
             --outfile=${outfile}             \
             --outformat="viewer"

echo
echo "######### written to ${outfile} ############"



