#-------------------------------------------------------------------------------------#
#                          visualize a structure with strands removed                 #
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
fn=simple
fn=fourhelix

# set the staples to retain.

if [ $fn == "simple" ]; then
   del_cmd="delete,retain=[5749504]"

elif [ $fn == "fourhelix" ]; then
    # staple colors are:
    #     1: [8,11184640],[41,243362
    #     0: [26,29184]  ,[60,7536862]
    #     3: [8, 5749504]
    #     2: [46, 243362 ], [60, 16225054] 
   #del_cmd="delete,retain=[11184640,243362]"
   del_cmd="delete,retain=[11184640]"
   del_cmd="delete,retain=[243362]"
   del_cmd="delete,retain=[29184]"
fi

# set the command to show all stand paths.
cmds="strand name=All  rep=path  show=true"

# execute the visualizer 
../../scripts/vis.py --infile=${dir}/${fn}.json       \
                     --inseqname=${seqname}           \
                     --atomic_model=${atomic_model}   \
                     --staples=${del_cmd}             \
                     --commands="${cmds}" 

