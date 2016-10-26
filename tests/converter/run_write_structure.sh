#-------------------------------------------------------------------------------------#
#                       write a dna structure to a json formated file                 #
#-------------------------------------------------------------------------------------#
# Set the cadnano design files directory.
data=../samples

fn=slottedcross
fn=Nature09_squarenut_no_joins
fn=square-test-1
fn=exampleOverhang
fn=fourhelix

outfile=./results/${fn}_structure.json

if [ ! -d "results/" ]; then
    mkdir results
fi


../../scripts/converter.py --infile=${data}/${fn}.json      \
             --inseqname=${seq}               \
             --informat="cadnano"             \
             --outfile=${outfile}             \
             --outformat="structure"


