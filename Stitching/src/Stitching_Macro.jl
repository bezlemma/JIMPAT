# This script runs multiple stitching jobs by directly calling functions from the Stitching module.
# Options include t_range z_range c_range save_location bin_factor align_strategy fusion_strategy feather_blend_size

include("Stitching.jl")
using .Stitching

# Example 1:
#job1_file = raw"C:\Users\uComp\Downloads\sample\20231113_1101_129_e14_kidney_w1_c_dapi_dba_ck8_sox9_10x2.companion.ome"
#stitch(loadfile=job1_file, bin_factor=2)

# Example 2: If you do not provide a load file, a GUI will pop up asking for a load file
#stitch()