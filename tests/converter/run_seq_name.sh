# This script tests converting a caDNAno design file into an Autodesk Nanodesign
# viewer file. 
#
# The sequence to assign to the design is given by a sequence name.

# Set the cadnano design files directory.
data=../samples

# Set the input caDNAno design file.
fn=fourhelix

# Set the sequence name. 
seq=M13mp18

# Set the sequence name.
# Valid names are: p7308, p7704, p8064, p8100, p8634, M13KO7, p7560 and M13mp18.
seq=M13mp18

# Modify the design with the inserts and deletes given in the caDNAno design file.
# Set to "true" to modify the design.
modify="false"

outfile=./results/${fn}_viewer.json 

if [ ! -d "results/" ]; then
    mkdir results
fi

../../scripts/converter.py --infile=${data}/${fn}.json  \
                           --informat="cadnano"         \
                           --inseqname=${seq}           \
                           --modify=${modify}           \
                           --outfile=${outfile}         \
                           --outformat="viewer"

