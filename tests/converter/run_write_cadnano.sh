# This script tests writing out a caDNAno design file.

# Set the cadnano design files directory.
data=../samples

# Set the input caDNAno design file.
fn=fourhelix

if [ ! -d "results/" ]; then
    mkdir results
fi

../../scripts/converter.py --infile=${data}/${fn}.json        \
                           --informat="cadnano"               \
                           --outfile=./results/${fn}_nd.json  \
                           --outformat="cadnano"

