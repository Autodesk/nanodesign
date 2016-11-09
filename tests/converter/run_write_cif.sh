# This script tests generating an atomistic structure from a caDNAno design file and
# writing it to a CIF-formatted file. 

# Set the cadnano design files directory.
data=../samples

# Set the input caDNAno design file.
fn=fourhelix

# Set the sequence name.
seq=M13mp18

outfile=./results/${fn}.cif

if [ ! -d "results/" ]; then
    mkdir results
fi

../../scripts/converter.py \
    --informat="cadnano"          \
    --infile=${data}/${fn}.json   \
    --inseqname=${seq}            \
    --outfile=${outfile}          \
    --outformat="cif"

