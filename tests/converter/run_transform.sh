#-------------------------------------------------------------------------------------#
#                     run the converter with a helix transformation.                  #
#-------------------------------------------------------------------------------------#
data=../samples

fn=fourhelix
transform="helices(0,1):rotate(90,0,0),translate(0,0,0)"
transform="helices(0,1):rotate(90,0,0),translate(0.5,0,0);helices(2,3):rotate(0,90,0),translate(0,0,0)"
transform="helices(1):rotate(0,0,90),translate(0,0,0)"
transform="helices(0,1):rotate(0,0,90),translate(0,0,0)"
transform="helices(0,1-3):rotate(0,0,90),translate(0,0,0)"

fn=slottedcross
transform="helices(0-41):connectors(scaffold)"

outfile=./results/${fn}_viewer.json  

if [ ! -d "results/" ]; then
    mkdir results
fi

../../scripts/converter.py --infile=${data}/${fn}.json \
             --informat="cadnano"        \
             --transform=${transform}    \
             --outfile=${outfile}        \
             --outformat="viewer"

