#-------------------------------------------------------------------------------------#
#                     write a dna stucture topology file                              #
#-------------------------------------------------------------------------------------#
# Set the cadnano design files directory.
data=../../tests/samples

fn=slottedcross
fn=square-test-1
fn=exampleOverhang
fn=Nature09_squarenut_no_joins
fn=fourhelix

converter.py --infile=${data}/${fn}.json              \
             --informat="cadnano"                     \
             --outfile=./results/${fn}_topology.json  \
             --outformat="topology"

