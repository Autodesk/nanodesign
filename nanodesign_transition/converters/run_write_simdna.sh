data=/Users/parkerda/nanodesign/test-data

fn=fourhelix
fn=slottedcross
fn=Rothemund-rect

out=/Users/parkerda/software/SimDNA/vis

converter.py --infile=${data}/${fn}.json --informat="cadnano"  --outfile=${out}/${fn}_nd.pairs  --outformat="simdna"

