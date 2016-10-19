data=../../tests/samples

fn=slottedcross
fn=Rothemund-rect
fn=Rothemund-rect_adjusted
fn=fourhelix-deletes
fn=fourhelix

outfile=./results/${fn}_nd.pairs
converter.py --infile=${data}/${fn}.json --informat="cadnano"  --modify=true --outfile=${outfile}  --outformat="simdna"

