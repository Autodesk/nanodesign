# This script tests converting a caDNAno design file into an Autodesk Nanodesign 
# viewer file. No sequence information is created.

# Set the cadnano design file directory.
data=../samples

# Set the input caDNAno design file.
fn=fourhelix

if [ ! -d "results/" ]; then
    mkdir results
fi

# Set the name of the output viewer file.
outfile=./results/${fn}_viewer.json  

../../scripts/converter.py --infile=${data}/${fn}.json \
                           --informat="cadnano"        \
                           --outfile=${outfile}        \
                           --outformat="viewer"

