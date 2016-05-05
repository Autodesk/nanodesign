
fn=slottedcross
fn=square-test-1
fn=exampleOverhang
fn=fourhelix
fn=Nature09_squarenut_no_joins

converter.py --infile=../data/${fn}.json --informat="cadnano"  --outfile=./results/${fn}_topology.json  --outformat="topology"

