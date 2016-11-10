# This shell script demonstrates how to execute the visualizer with command-line commands.

# Set the directory path to cadnano design files.
dir=../samples/

# Set the sequence name.
seqname=M13mp18

# Set the cadnano file name.
fn=fourhelix

# Set the flag for generating atomic structures.
atomic_model="true"

# Examples of some visualizer commands.
if [ $fn == "fourhelix" ]; then
   cmds="helix name=0 rep=geometry color=(1,0.5,0) show=true"
   cmds="strand name=staple_0_26  rep=path  color=(1,0,0) show=true"
   cmds="strand names=start_helices[0]  rep=path  color=(1,0,0) line_width=4.0  show=true"
   cmds="strand names=in_helices[19]  rep=path  color=(1,0,0) line_width=4.0  show=true"
fi

../../scripts/vis.py --infile=${dir}/${fn}.json     \
                     --inseqname=${seqname}         \
                     --atomic_model=${atomic_model} \
                     --commands="${cmds}"

