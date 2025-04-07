module Stitching

using FileIO,Printf, FFTW, Base.Threads, LinearAlgebra, TiffImages
using Statistics, SparseArrays, OrderedCollections
using NativeFileDialog, ProgressMeter
using EzXML, DataFrames

include("Stitching_Load_Functions.jl")
include("Stitching_Fuse_Functions.jl")
include("Stitching_Align_Functions.jl")
include("Stitching_GUI.jl")
include("SaveTiffPlz.jl")

function main(companion_file::String, t_range::AbstractRange, z_range::AbstractRange, c_range::AbstractRange)

    # --- Parameters ---
    align_strategy = "From Position" #TODO: Implement global alginment
    bin_factor = 4
    Fusion_Strategy = "None" #TODO: Re-do "Feather"
    feather_blend_size = 30.0


    println("Using Companion File: ", companion_file)
    base_path = dirname(companion_file)
    image_data = parse_ome(companion_file)
    meta_data = parse_ome_metadata(companion_file)
    # --- Validate User-Provided Ranges ---
    metadata_size_t = meta_data["SizeT"]
    metadata_size_z = meta_data["SizeZ"]
    metadata_size_c = meta_data["SizeC"]

    if maximum(t_range) > metadata_size_t || minimum(t_range) < 1
        println("Error: Provided Time range ($(t_range)) is outside the valid metadata range (1:$(metadata_size_t)).")
        return
    end
    if maximum(z_range) > metadata_size_z || minimum(z_range) < 1
        println("Error: Provided Z range ($(z_range)) is outside the valid metadata range (1:$(metadata_size_z)).")
        return
    end
    if maximum(c_range) > metadata_size_c || minimum(c_range) < 1
        println("Error: Provided Channel range ($(c_range)) is outside the valid metadata range (1:$(metadata_size_c)).")
        return
    end

    #Alignment
    if align_strategy == "From Position"
        offsets = tile_using_positions(image_data, bin_factor)
    else
        println("Invalid alignment stratgy, exiting.")
        return
    end

    #Fusion
    println("\nStarting fusion process...")
    stitched_series = fuse_main(image_data, offsets, Fusion_Strategy, base_path,
        c_range, z_range, t_range, bin_factor,
        feather_blend_size)

    # --- Saving Code  ---
    binned_str = bin_factor > 1 ? "_binned$(bin_factor)x" : ""
    fusion_str = Fusion_Strategy != "None" ? "_$(lowercase(Fusion_Strategy))" : ""
    output_filename = replace(basename(companion_file), r"(\.companion)?\.ome$"i => "") * "_stitched$(binned_str)$(fusion_str).tif"
    output_file = joinpath(base_path, output_filename)
    println("Saving stitched 5D array to: $output_file")
    stitched_16 = clamp.(stitched_series, 0.0f0, 1.0f0) .* 65535.0f0
    stitched_16 = round.(UInt16, stitched_16)
    SaveTiffPlz(output_file, stitched_16)
end

function main(companion_file::String)
    meta_data = parse_ome_metadata(companion_file)
    t_range = 1:meta_data["SizeT"]
    z_range = 1:meta_data["SizeZ"]
    c_range = 1:meta_data["SizeC"]
    main(companion_file, t_range, z_range, c_range)
end

function julia_main()::Cint
    println("Opening file selection dialog...")
    filter_list = "companion.ome"
    selected_file = ""
    selected_file = pick_file(pwd(); filterlist=filter_list)
    if isempty(selected_file)
        return 1
    end
    main(selected_file)
    return 0
end

#Example for running with a given file and given ranges 
#      specific_t_range = 1:1
#      specific_z_range = 5:10  # Example: subset of Z
#      specific_c_range = 1:2   # Example: subset of C
#      main(companion_file_ex3, specific_t_range, specific_z_range, specific_c_range)

#Examples for running with a given file, all range
#base_path = raw"C:\Users\uComp\Documents\LinneaData\20250226\Sample5_1\\"
#companion_file = base_path * "20250226_ll_sample5_60x_10msexp_1-2bint_gfp-bypass_60sint1.companion.ome"

#base_path = raw"C:\Users\uComp\Downloads\sample\\"
#companion_file = base_path * "20231113_1101_129_e14_kidney_w1_c_dapi_dba_ck8_sox9_10x2.companion.ome"
#main(companion_file)


#base_path = raw"C:\Users\Bez\Downloads\stitching\data\sample\\"
#companion_file = base_path * "20231113_1101_129_e14_kidney_w1_c_dapi_dba_ck8_sox9_10x2.companion.ome"
#main(companion_file)


## Example for GUI run
#main()

end
#Debugging:
#  c_range = 1:1
#  z_range = 50:100
#  t_range = 1:1
#  bin_factor = 2
#  Fusion_Strategy = "None"
#  feather_blend_size = 30.0

 # image_data = parse_ome(companion_file)
#meta_data = parse_ome_metadata(companion_file)
#  offsets = tile_using_positions(image_data, bin_factor)
#  stitched_series = fuse_main(image_data, offsets, Fusion_Strategy, base_path,
#      c_range, z_range, t_range, bin_factor,
#      feather_blend_size)