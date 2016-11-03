#-------------------------------------------------------------------------------------#
#                   visualize a structure with added maximal strands                  #
#-------------------------------------------------------------------------------------#

# directory path to cadnano design files.
dir=../samples/

# sequence file name
seqfile=fourhelix_M13mp18

# sequence name
seqname=M13mp18

# set the flag for generating atomic structures.
atomic_model="true"
atomic_model="false"

# set cadnano file name.
fn=fourhelix

# set the staples to retain.

if [ $fn == "simple" ]; then
   strand_cmd="maximal_set,retain=[5749504]"

elif [ $fn == "fourhelix" ]; then
    # staple colors are:
    #     1: [8,11184640],[41,243362
    #     0: [26,29184]  ,[60,7536862]
    #     3: [8, 5749504]
    #     2: [46, 243362 ], [60, 16225054] 
   strand_cmd="maximal_set,retain=[11184640,243362]"
   strand_cmd="maximal_set,retain=[29184]"
   strand_cmd="maximal_set"

    # set the command to show all stand paths.
    cmds="strand name=All  rep=path  show=true;helix name=0  rep=maximal_crossovers  show=true"
fi

# execute the visualizer 
../../scripts/vis.py --infile=${dir}/${fn}.json       \
                     --inseqname=${seqname}           \
                     --atomic_model=${atomic_model}   \
                     --staples=${strand_cmd}          \
                     --commands="${cmds}" 

