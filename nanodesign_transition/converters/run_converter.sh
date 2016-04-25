
fn=slottedcross
fn=fourhelix
fn=Nature09_squarenut_no_joins
fn=square-test-1
fn=exampleOverhang

converter.py --infile=../data/${fn}.json --informat="cadnano"  --outfile=./results/${fn}_viewer.json  --outformat="viewer"

