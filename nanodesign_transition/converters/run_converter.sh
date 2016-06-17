
data=../../tests/samples

fn=slottedcross
fn=Nature09_squarenut_no_joins
fn=square-test-1
fn=exampleOverhang
fn=hc-test-9
fn=protractor_180_98_v3
fn=fourhelix
fn=42hb

outfile=./results/${fn}_viewer.json  

converter.py --infile=${data}/${fn}.json \
             --informat="cadnano"        \
             --outfile=${outfile}        \
             --outformat="viewer"

