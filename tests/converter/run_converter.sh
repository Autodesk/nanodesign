data=../samples

fn=42hb
fn=exampleOverhang
fn=fourhelix
fn=hc-test-9
fn=nanorobot.v2
fn=Nature09_squarenut_no_joins
fn=protractor_180_98_v3
fn=slottedcross
fn=square-test-1

fn=fourhelix

if [ ! -d "results/" ]; then
    mkdir results
fi

outfile=./results/${fn}_viewer.json  

../../scripts/converter.py --infile=${data}/${fn}.json \
                           --informat="cadnano"        \
                           --outfile=${outfile}        \
                           --outformat="viewer"

