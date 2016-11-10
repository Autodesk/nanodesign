# This shell script demonstrates how to execute the visualizer to show the maximal staple set
# for a design.
#
# A list of staple colors, matching those defined for staples in the cadnano files, may
# be given to identify the staples retained after the maximal staple set for a design
# is calculated.

# Set the directory path to cadnano design files.
dir=../samples/

# Set the sequence name.
seqname=M13mp18

# Set the flag for generating atomic structures.
atomic_model="false"

# Set cadnano file name.
fn=fourhelix

# Set the strand operation command.

if [ $fn == "fourhelix" ]; then
    # The valid staple colors from the fourhelix design file (helix number: [position,color], ...):
    #     1: [8,11184640],[41,243362
    #     0: [26,29184]  ,[60,7536862]
    #     3: [8, 5749504]
    #     2: [46, 243362 ], [60, 16225054] 
    strand_cmd="maximal_set,retain=[11184640,243362]"
    strand_cmd="maximal_set,retain=[29184]"
    strand_cmd="maximal_set"
    # Set the command to show all stand paths.
    cmds="strand name=All  rep=path  show=true;helix name=0  rep=maximal_crossovers  show=true"
else
    strand_cmd="maximal_set"
    cmds="strand name=All  rep=path  show=true"
fi

../../scripts/vis.py --infile=${dir}/${fn}.json       \
                     --inseqname=${seqname}           \
                     --atomic_model=${atomic_model}   \
                     --staples=${strand_cmd}          \
                     --commands="${cmds}" 
