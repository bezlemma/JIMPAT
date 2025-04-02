using FileIO, EzXML, FFTW, Base.Threads, LinearAlgebra, TiffImages
using Statistics, SparseArrays, OrderedCollections
using NativeFileDialog, ProgressMeter

include("Stitching_Load_Functions.jl")
include("Stitching_Fuse_Functions.jl")
include("Stitching_Align_Functions.jl")
include("SaveTiffPlz.jl")

# --- Main Execution Logic ---

function main(companion_file::String)
    println("Using provided companion file: ", companion_file)
    base_path = dirname(companion_file)
    println("Parsing OME companion file...")
    image_data = parse_ome_companion(companion_file)
    if isempty(image_data)
        println("Error: Parsing failed.")
        return
    end
    println("Parsed data for $(length(image_data)) valid tiles.")
    if !haskey(image_data[1], "metadata_size_c")
        println("Error: Cannot determine metadata dimensions.")
        return
    end
    metadata_size_t = image_data[1]["metadata_size_t"]
    metadata_size_z = image_data[1]["metadata_size_z"]
    metadata_size_c = image_data[1]["metadata_size_c"]
    println("Dimensions from Metadata: T=$metadata_size_t, Z=$metadata_size_z, C=$metadata_size_c")
    if metadata_size_t <= 0 || metadata_size_z <= 0 || metadata_size_c <= 0
        println("Error: Metadata dimensions invalid.")
        return
    end

    align_strategy = "From Position"
    time_range = 1:metadata_size_t
    z_range = 1:metadata_size_z
    c_range = 1:metadata_size_c
    bin_factor = 1
    Fusion_Strategy = "None"
    feather_blend_size = 30.0

    println("Calculating tile offsets...")
    offsets = tile_using_positions(image_data, bin_factor)
    if isempty(offsets)
        println("Error: No offsets calculated.")
        return
    end
    println("Calculated offsets for $(length(offsets)) tiles.")

    stitched_series = fuse_main(image_data, offsets, Fusion_Strategy, base_path, c_range, z_range, time_range, bin_factor, metadata_size_c, metadata_size_z, metadata_size_t, feather_blend_size)
    if isempty(stitched_series) || all(size(stitched_series) .== 0)
        println("Error: Stitching returned empty result.")
        return
    end

    #Saving Code
    binned_str = bin_factor > 1 ? "_binned$(bin_factor)x" : ""
    fusion_str = Fusion_Strategy != "None" ? "_$(lowercase(Fusion_Strategy))" : ""
    output_filename = replace(basename(companion_file), r"(\.companion)?\.ome$"i => "") * "_stitched$(binned_str)$(fusion_str).tif"
    output_file = joinpath(base_path, output_filename)
    println("Saving stitched 5D array to: $output_file")
    stitched_16 = clamp.(stitched_series, 0.0f0, 1.0f0) .* 65535.0f0
    stitched_16 = round.(UInt16, stitched_16)
    TiffSaver.SaveTiffPlz(output_file, stitched_16)
    
    println("\nDone.")
end

# Version of main with no arguments - triggers file dialog
function main()
    println("No companion file path provided."); println("Opening file selection dialog...")
    filter_list = "companion.ome"
    try; selected_file = pick_file(pwd(); filterlist=filter_list)
        if isempty(selected_file); println("File selection cancelled. Exiting."); return; end
        println("File selected: ", selected_file); main(selected_file)
    catch e; println("\nError during file selection dialog:"); end
end

## Code run examples
main()

#base_path = raw"C:\Users\uComp\Documents\LinneaData\20250226\Sample5_1\\"
#companion_file = base_path * "20250226_ll_sample5_60x_10msexp_1-2bint_gfp-bypass_60sint1.companion.ome"

#base_path = raw"C:\Users\uComp\Downloads\sample\\"
#companion_file = base_path * "20231113_1101_129_e14_kidney_w1_c_dapi_dba_ck8_sox9_10x2.companion.ome"
#main(companion_file)

