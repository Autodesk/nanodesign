
dir=/Users/parkerda/nanodesign/parkerda/nanodesign/nanodesign_transition/converters

${dir}/converter.py --infile=${1}.json --informat="cadnano"  --outfile=${1}_viewer.json     --outformat="viewer"

${dir}/converter.py --infile=${1}.json --informat="cadnano"  --outfile=${1}_topology.json   --outformat="topology"

${dir}/converter.py --infile=${1}.json --informat="cadnano"  --outfile=${1}_structure.json  --outformat="structure"



