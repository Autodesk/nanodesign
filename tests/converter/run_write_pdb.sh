
data=../samples

seq=M13mp18

fn=42hb
fn=AutodeskA
fn=doublegear
fn=endcap_issues
fn=exampleOverhang
fn=fourhelix
fn=fourhelix-1
fn=fourhelix-2
fn=hc-test
fn=hc-test-1
fn=hc-test-3
fn=hc-test-6    # does not show in chimera
fn=hc-test-8    # process single stranded, does not show in chimera.
fn=hc-test-10
fn=monolith
fn=Nature09_squarenut_no_joins
fn=nanorobot.v2
fn=pointer
fn=robot_v1.9_bent_2
fn=simple
fn=simple-1
fn=six_helix_honeycomb
fn=slottedcross
fn=sq-test-1
fn=sq-test-2
fn=sq-test-3
fn=sq-test-4
fn=sq-test-5
fn=sq-test-6
fn=square-test-1
fn=square-test-2
fn=square-test-3

fn=fourhelix
fn=Rothemund-rect

outfile=./results/${fn}.pdb

if [ ! -d "results/" ]; then
    mkdir results
fi

../../scripts/converter.py \
    --informat="cadnano"          \
    --infile=${data}/${fn}.json   \
    --inseqname=${seq}            \
    --outfile=${outfile}          \
    --outformat="pdb"

#cp ${outfile} ~ 

