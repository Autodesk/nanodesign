#-------------------------------------------------------------------------------------#
#                                 run the auto-stapler                                #
#-------------------------------------------------------------------------------------#
# Set the cadnano design files directory.
data=../samples

# Execute the auto-stapler.
fn=6hb
fn=42hb
fn=rr_triangle

python ../../scripts/stapler.py ${data}/${fn}.json  


