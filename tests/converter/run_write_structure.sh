# This script tests writing a dna structure created from a caDNAno design file to JSON-formated 
# and plain text files. The files contain information about the bases, strands and domains created
# from a design.

# Set the cadnano design files directory.
data=../samples

# Set the input caDNAno design file.
fn=fourhelix

outfile=./results/${fn}_structure.json

if [ ! -d "results/" ]; then
    mkdir results
fi


../../scripts/converter.py --infile=${data}/${fn}.json  \
                           --inseqname=${seq}           \
                           --informat="cadnano"         \
                           --outfile=${outfile}         \
                           --outformat="structure"
