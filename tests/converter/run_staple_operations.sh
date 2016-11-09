#-------------------------------------------------------------------------------------#
#                         run the converter with staple operations                    #
#-------------------------------------------------------------------------------------#
data=../samples

# Set caDNAno design file.
fn=fourhelix

# Set operation.
op=delete
op=maximal_set

# Set commands depending on input file name.
if [ $fn == "fourhelix" ]; then
    # Delete staples.
    if [ $op == "delete" ]; then
        del_cmd="delete,retain=[243362]"
        del_cmd="delete"
        del_cmd="delete,retain=[243362,29184,7536862,5749504,16225054]"
        del_cmd="delete,retain=[11184640]"
        del_cmd="delete,retain=[29184]"
        del_cmd="delete,retain=[11184640,243362]"
        staple_cmd=${del_cmd}
        outfile=./results/${fn}_nd_del.json  
    fi

    # Create maximal staple set.
    if [ $op == "maximal_set" ]; then
        mset_cmd="maximal_set,retain=[11184640]"
        mset_cmd="maximal_set,retain=[243362]"
        mset_cmd="maximal_set,retain=[29184]"
        mset_cmd="maximal_set,retain=[11184640,243362]"
        mset_cmd="maximal_set"
        staple_cmd=${mset_cmd}
        outfile=./results/${fn}_nd_max.json  
    fi
fi

if [ ! -d "results/" ]; then
    mkdir results
fi

../../scripts/converter.py --infile=${data}/${fn}.json   \
                           --informat="cadnano"          \
                           --staples=${staple_cmd}       \
                           --outfile=${outfile}          \
                           --outformat="cadnano"

