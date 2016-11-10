# This script tests converting a caDNAno design file into an Autodesk Nanodesign
# viewer file. 
#
# The sequence to assign to the design is given by a CSV file with the following
# naming convention:
#
#    <caDNAnoFileName>_<SequenceName>.csv 

# Set the cadnano design files directory.
data=../samples

# Set the input caDNAno design file.
fn=fourhelix

# Set the sequence name. 
seq=M13mp18

outfile=./results/${fn}_viewer.json 

if [ ! -d "results/" ]; then
    mkdir results
fi

../../scripts/converter.py --infile=${data}/${fn}.json             \
             --inseqfile=${data}/${fn}_${seq}.csv    \
             --informat="cadnano"                    \
             --outfile=${outfile}                    \
             --outformat="viewer"

