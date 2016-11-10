# This script tests writing a dna structure created from a caDNAno design file to JSON-formated 
# and plain text files. The files contain information about the bases created from a design.

# Set the cadnano design files directory.
data=../samples

# Set the input caDNAno design file.
fn=fourhelix

if [ ! -d "results/" ]; then
    mkdir results
fi

../../scripts/converter.py --infile=${data}/${fn}.json              \
                           --informat="cadnano"                     \
                           --outfile=./results/${fn}_topology.json  \
                           --outformat="topology"

