#-------------------------------------------------------------------------------------#
#                     run the converter with a sequence name.                         #
#-------------------------------------------------------------------------------------#

# Set the cadnano design files directory.
data=../../tests/samples

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
fn=fourhelix_inserts
fn=fourhelix_deletes
fn=fourhelix_inserts_2
fn=fourhelix_ss_inserts
fn=protractor_30_98_v4
fn=fourhelix
fn=slottedcross

outfile=./results/${fn}_viewer.json 
#outfile=/Users/parkerda/tirrenu/res/testDNA/fourhelix_vis.json 


converter.py --infile=${data}/${fn}.json      \
             --informat="cadnano"             \
             --inseqname=${seq}               \
             --modify=${modify}               \
             --outfile=${outfile}             \
             --outformat="viewer"

echo
echo "######### written to ${outfile} ############"

compare_script=/Users/parkerda/nanodesign/parkerda-latest/nanodesign/scripts/compare-viewer-json.py 

gold=/Users/parkerda/software/nanodesign/testing/viewer/results/gold/
gold=/Users/parkerda/nanodesign/parkerda-gold/nanodesign/nanodesign/nanodesign_transition/converters/results/

${compare_script} ${gold}/${fn}_viewer.json ${outfile} &> out

if (($? == 0)); then
   echo "Comparison passed. "
else
   echo "**** Comparison failed." 
fi


