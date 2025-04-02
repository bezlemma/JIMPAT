using FileIO, EzXML, FFTW, Base.Threads, LinearAlgebra, TiffImages
using Statistics, SparseArrays, OrderedCollections
using NativeFileDialog, ProgressMeter
using GLMakie

# Assuming these files contain the necessary function definitions
include("Stitching_Load_Functions.jl")
include("Stitching_Fuse_Functions.jl")
include("Stitching_Align_Functions.jl")
include("Stitching_GUI.jl")
include("SaveTiffPlz.jl")

function main(companion_file::String, time_range::AbstractRange, z_range::AbstractRange, c_range::AbstractRange)
    println("Using Companion File: ", companion_file)
    base_path = dirname(companion_file)
    image_data = parse_ome_companion(companion_file)
    println("Parsed data for $(length(image_data)) valid tiles.")
    metadata_size_t = image_data[1]["metadata_size_t"]
    metadata_size_z = image_data[1]["metadata_size_z"]
    metadata_size_c = image_data[1]["metadata_size_c"]

    # --- Validate User-Provided Ranges ---
    if maximum(time_range) > metadata_size_t || minimum(time_range) < 1
        println("Error: Provided Time range ($(time_range)) is outside the valid metadata range (1:$(metadata_size_t)).")
        return;    end
    if maximum(z_range) > metadata_size_z || minimum(z_range) < 1
        println("Error: Provided Z range ($(z_range)) is outside the valid metadata range (1:$(metadata_size_z)).")
        return;    end
    if maximum(c_range) > metadata_size_c || minimum(c_range) < 1
        println("Error: Provided Channel range ($(c_range)) is outside the valid metadata range (1:$(metadata_size_c)).")
        return;    end

    # --- Parameters ---
    align_strategy = "From Position"
    bin_factor = 2
    Fusion_Strategy = "None"
    feather_blend_size = 30.0

    offsets = tile_using_positions(image_data, bin_factor)

    println("\nStarting fusion process...")
    stitched_series = fuse_main(image_data, offsets, Fusion_Strategy, base_path,
        c_range, z_range, time_range, bin_factor,
        metadata_size_c, metadata_size_z, metadata_size_t,
        feather_blend_size)

    # --- Saving Code  ---
    binned_str = bin_factor > 1 ? "_binned$(bin_factor)x" : ""
    fusion_str = Fusion_Strategy != "None" ? "_$(lowercase(Fusion_Strategy))" : ""
    output_filename = replace(basename(companion_file), r"(\.companion)?\.ome$"i => "") * "_stitched$(binned_str)$(fusion_str).tif"
    output_file = joinpath(base_path, output_filename)
    println("Saving stitched 5D array to: $output_file")
    stitched_16 = clamp.(stitched_series, 0.0f0, 1.0f0) .* 65535.0f0
    stitched_16 = round.(UInt16, stitched_16)
    TiffSaver.SaveTiffPlz(output_file, stitched_16)
end

function main(companion_file::String)
    println("Using Companion File: ", companion_file)
    image_data = parse_ome_companion(companion_file)
    metadata_size_t = image_data[1]["metadata_size_t"]
    metadata_size_z = image_data[1]["metadata_size_z"]
    metadata_size_c = image_data[1]["metadata_size_c"]
    main(companion_file, 1:metadata_size_t, 1:metadata_size_z, 1:metadata_size_c)
end

function main()
    println("Opening file selection dialog...")
    filter_list = "companion.ome"
    selected_file = ""
    selected_file = pick_file(pwd(); filterlist=filter_list)
    if isempty(selected_file)
        println("File selection cancelled. Exiting.")
        return
    end
    println("File selected: ", selected_file)

    local image_data
    image_data = parse_ome_companion(selected_file)
    metadata_size_t = image_data[1]["metadata_size_t"]
    metadata_size_z = image_data[1]["metadata_size_z"]
    metadata_size_c = image_data[1]["metadata_size_c"]

    println("\nOpening GUI to select ranges (using Makie)...")
    gui_result = prompt_for_ranges_makie(metadata_size_t, metadata_size_z, metadata_size_c)

    if isnothing(gui_result)
        println("Range selection cancelled or GUI closed. Exiting.")
        return
    else
        time_range, z_range, c_range = gui_result
        println("\nProceeding with selected file and ranges from GUI...")
        main(selected_file, time_range, z_range, c_range) # Call the core function
    end
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


## Example for GUI run
main()



