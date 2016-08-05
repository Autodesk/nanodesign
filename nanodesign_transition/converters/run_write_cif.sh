
data=../../tests/samples

fn=square-test-1
fn=exampleOverhang
fn=endcap_issues
fn=simple
fn=Rothemund-rect
fn=monolith
fn=doublegear
fn=pointer
fn=sq-test-1
fn=Nature09_squarenut_no_joins
fn=slottedcross
fn=hc-test
fn=aNANO_3D_7_14_final
fn=fourhelix
fn=robot_v1.9_bent_2

fn=nanorobot.v2


seq=M13mp18

outfile=./results/${fn}.cif

converter.py \
    --informat="cadnano"          \
    --infile=${data}/${fn}.json   \
    --inseqname=${seq}            \
    --outfile=${outfile}          \
    --outformat="cif"

cp ${outfile} ~ 

