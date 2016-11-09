#-------------------------------------------------------------------------------------#
#                          visualize a transformed structure                          #
#-------------------------------------------------------------------------------------#

# directory path to cadnano design files.
dir=../samples/

# sequence file name
seqfile=fourhelix_M13mp18

# sequence name
seqname=M13mp18

# set cadnano file name
fn=convex_triangle

# set the flag for generating atomic structures.
atomic_model="true"
atomic_model="false"

if [ $fn == "slottedcross" ]; then 
    transform="helices(0-41):connectors(scaffold)"
    transform="helices(0-41):rotate(270,0,0),translate(0,0,41.4)"
    cmds="strand name=Scaffold_55_7 rep=path show=true;strand name=Scaffold_55_7 rep=connectors show=true;helix name=19  rep=geometry show=true;helix name=46 rep=geometry show=true;helix name=19  rep=base_positions  show=true"

elif [ $fn == "convex_triangle" ]; then
    transform="helices(1-10):connectors(scaffold)"
   transform=""
else
   transform=""
fi

# execute the visualizer 
../../scripts/vis.py --infile=${dir}/${fn}.json       \
                     --inseqname=${seqname}           \
                     --atomic_model=${atomic_model}   \
                     --transform=${transform}         \
                     --commands="${cmds}" 

