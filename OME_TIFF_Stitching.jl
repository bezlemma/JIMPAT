using FileIO, FFTW, Base.Threads, LinearAlgebra, TiffImages
using Statistics
using SparseArrays  # Needed for the global solver
using OrderedCollections

include("./Stitching_Load_Functions.jl")
include("./Stitching_Align_Functions.jl")
include("./Stitching_Blend_Functions.jl")

#Load the data.
base_path = raw"C:\Users\uComp\Documents\LinneaData\20250226\Sample5_1\\"
companion_file = base_path * "20250226_ll_sample5_60x_10msexp_1-2bint_gfp-bypass_60sint1.companion.ome"

#base_path = raw"C:\Users\uComp\Downloads\sample\\"
#companion_file = base_path * "20231113_1101_129_e14_kidney_w1_c_dapi_dba_ck8_sox9_10x2.companion.ome"

image_data = parse_ome_companion(companion_file)
size_t = image_data[1]["size_t"]; size_z = image_data[1]["size_z"]; size_c = image_data[1]["size_c"];

## Align the images
#TODO: Add Option to align off a merge of AmodB frames
#TODO: Add Option to do a new alignment for every frame

align_strategy = "From Position"; #"From Position, From Frame N, From Nth Frames, For Each Frame
Frame_N = 42; #TODO: Right now this is time, but we need it to be time or color or z
time_range = 1:size_t #40:45
z_range = 1:size_z 
c_range = 1:size_c
bin_factor = 2; #Bin1,2 or 4 supported. Defaults to 2
background_subtract = false; #TODO: We should allow the option of background subtraction here if the user wants.

if align_strategy == "From Position"
    offsets = tile_using_positions(image_data, bin_factor)
elseif align_strategy == "From Frame N"
    # if Frame_N > (size_t - 1); Frame_N = (size_t - 1); end
    # println("Aligning tiles based on final timepoint (index $Frame_N).")

    # #TODO This should a global alignment
    # (offsets, rowcol_image) = align_tiles_reference_timepoint(image_data, base_path;
    #     time_index=Frame_N, overlap_frac=0.1, max_allowed_shift=200.0) 

    # println("Done computing offsets.")
end

## 3 is to blend the images
Fusion_Strategy = "None" #"None", "Feather"
stitched_series = stitch_main(image_data, offsets, Fusion_Strategy, base_path,time_range,z_range,c_range)

## Save the file
println("Saving...")
output_file = raw"C:\Users\uComp\Documents\LinneaData\20250226\Sample5_1\output_Julia_binned.tif"
stitched_16 = clamp.(stitched_series, 0f0, 1f0) .* 65535f0
stitched_16 = round.(UInt16, stitched_16)
save(output_file, stitched_16)
println("Saved stitched binned series (Float32) to: $output_file")
