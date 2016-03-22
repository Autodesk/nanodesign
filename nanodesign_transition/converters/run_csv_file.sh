#-------------------------------------------------------------------------------------#
#                     run the converter with a sequence file name.                    #
#-------------------------------------------------------------------------------------#

fn=fourhelix_skip
fn=fourhelix

seq=p7308
seq=p7560
seq=p7704
seq=p8064
seq=p8100
seq=p8634
seq=M13KO7
seq=M13mp18

converter.py --infile=../data/${fn}.json                  \
             --inseqfile=../data/${fn}_${seq}.csv         \
             --informat="cadnano"                         \
             --outfile=./results/${fn}_viewer.json        \
             --outformat="viewer"

