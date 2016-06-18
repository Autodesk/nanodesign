data=../../tests/samples

fn=slottedcross
fn=Rothemund-rect
fn=fourhelix

outfile=./results/${fn}_nd.pairs

converter.py --infile=${data}/${fn}.json --informat="cadnano"  --outfile=${outfile}  --outformat="simdna"

