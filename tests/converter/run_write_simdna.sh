# This script tests converting a caDNAno design file into a SimDNA file. 

# Set the cadnano design files directory.
data=../samples

# Set the input caDNAno design file.
fn=fourhelix

outfile=./results/${fn}_nd.pairs

if [ ! -d "results/" ]; then
    mkdir results
fi

../../scripts/converter.py --infile=${data}/${fn}.json   \
                           --informat="cadnano"          \
                           --modify=true                 \
                           --outfile=${outfile}          \
                           --outformat="simdna"
