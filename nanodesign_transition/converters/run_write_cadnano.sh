#-------------------------------------------------------------------------------------#
#                             write out a cadnano file                                #
#-------------------------------------------------------------------------------------#
data=../../tests/samples

fn=Nature09_squarenut_no_joins
fn=square-test-1
fn=exampleOverhang
fn=slottedcross
fn=Rothemund-rect
fn=robot_v1.9_bent_2
fn=endcap_issues
fn=fourhelix

converter.py --infile=${data}/${fn}.json        \
             --informat="cadnano"               \
             --outfile=./results/${fn}_nd.json  \
             --outformat="cadnano"

