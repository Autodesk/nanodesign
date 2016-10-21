#-------------------------------------------------------------------------------------#
#                       write a dna structure to a json formated file                 #
#-------------------------------------------------------------------------------------#
# Set the cadnano design files directory.
data=../../tests/samples

fn=slottedcross
fn=Nature09_squarenut_no_joins
fn=square-test-1
fn=exampleOverhang
fn=fourhelix

outfile=./results/${fn}_structure.json

converter.py --infile=${data}/${fn}.json      \
             --inseqname=${seq}               \
             --informat="cadnano"             \
             --outfile=${outfile}             \
             --outformat="structure"


