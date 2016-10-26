#-------------------------------------------------------------------------------------#
#                     write a dna stucture topology file                              #
#-------------------------------------------------------------------------------------#
# Set the cadnano design files directory.

data=../samples

fn=slottedcross
fn=square-test-1
fn=exampleOverhang
fn=Nature09_squarenut_no_joins
fn=fourhelix

if [ ! -d "results/" ]; then
    mkdir results
fi

../../scripts/converter.py --infile=${data}/${fn}.json              \
             --informat="cadnano"                     \
             --outfile=./results/${fn}_topology.json  \
             --outformat="topology"

