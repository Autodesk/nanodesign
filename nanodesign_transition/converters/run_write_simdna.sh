data=../../tests/samples

fn=slottedcross
fn=Rothemund-rect_adjusted
fn=fourhelix-deletes
fn=fourhelix-deletes-3
fn=Rothemund-rect
fn=protractor_30_98_v4
fn=protractor_30_no_wrap
fn=Rothemund-rect_adjusted
fn=monolith_right_twist_no_wrap
fn=monolith_right_twist_no_wrap_no_ss
fn=fourhelix

outfile=./results/${fn}_nd.pairs

converter.py --infile=${data}/${fn}.json   \
             --informat="cadnano"          \
             --modify=true                 \
             --outfile=${outfile}          \
             --outformat="simdna"

