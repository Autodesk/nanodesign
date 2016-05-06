
fn=slottedcross
fn=Nature09_squarenut_no_joins
fn=square-test-1
fn=exampleOverhang
fn=fourhelix

converter.py --infile=../data/${fn}.json --informat="cadnano"  --outfile=./results/${fn}_structure.json  --outformat="structure"

