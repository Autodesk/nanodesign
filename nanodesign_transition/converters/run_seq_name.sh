#-------------------------------------------------------------------------------------#
#                     run the converter with a sequence name.                         #
#-------------------------------------------------------------------------------------#
data=../../tests/samples

fn=fourhelix_skip
fn=protractor_180_98_v3
fn=robot_v1.9_bent_2
fn=fourhelix
fn=42hb

seq=p7308
seq=p7704
seq=p8064
seq=p8100
seq=p8634
seq=M13KO7
seq=p7560
seq=M13mp18

outfile=./results/${fn}_viewer.json 

converter.py --infile=${data}/${fn}.json      \
             --inseqname=${seq}               \
             --informat="cadnano"             \
             --outfile=${outfile}             \
             --outformat="viewer"

