using FileIO, FFTW, Base.Threads, LinearAlgebra, TiffImages
using Statistics
using SparseArrays
using OrderedCollections
using NativeFileDialog # <-- Add this import

include("./Stitching_Load_Functions.jl")
include("./Stitching_Align_Functions.jl")
include("./Stitching_Blend_Functions.jl")

# Version of main that takes a file path string
function main(companion_file::String)
    println("Using provided companion file: ", companion_file)

    # --- Start of your original main function logic ---
    image_data = parse_ome_companion(companion_file)
    if isempty(image_data) # Add a check in case parsing fails or file is invalid
        println("Error: Could not parse OME companion file or it's empty.")
        return
    end
    size_t = image_data[1]["size_t"]
    size_z = image_data[1]["size_z"]
    size_c = image_data[1]["size_c"]

    ## Align the images
    align_strategy = "From Position" #"From Position, From Frame N, From Nth Frames, For Each Frame
    Frame_N = 42 # Example value, adjust as needed
    time_range = 1:size_t
    z_range = 1:size_z
    c_range = 1:size_c
    bin_factor = 2 # Bin 1, 2 or 4 supported. Defaults to 2
    background_subtract = false

    local offsets # Define offsets variable in the function scope

    if align_strategy == "From Position"
        offsets = tile_using_positions(image_data, bin_factor)
    elseif align_strategy == "From Frame N"
        # Placeholder for alignment logic using Frame_N
        println("Aligning tiles based on Frame N (logic needs implementation).")
        # Example: You would compute offsets here based on Frame_N
        # offsets = compute_offsets_for_frame_n(...)
        println("Warning: 'From Frame N' alignment strategy is not fully implemented.")
        # For now, let's fall back to position-based or handle error
        println("Falling back to position-based alignment.")
        offsets = tile_using_positions(image_data, bin_factor) # Example fallback
    else
         println("Error: Unknown alignment strategy: ", align_strategy)
         return # Exit if strategy is unknown
    end

    if !isdefined(Main, :offsets) || isempty(offsets) # Check if offsets were successfully computed
         println("Error: Failed to compute tile offsets.")
         return
    end

    ## Blend the images
    Fusion_Strategy = "None" #"None", "Feather"
    # Determine base_path from the companion_file path
    base_path = dirname(companion_file)
    stitched_series = stitch_main(image_data, offsets, Fusion_Strategy, base_path, time_range, z_range, c_range)

    ## Save the file
    println("Saving...")
    # Construct output path relative to the input file's directory
    output_filename = replace(basename(companion_file), ".companion.ome" => "_stitched_binned.tif")
    output_file = joinpath(base_path, output_filename)

    stitched_16 = clamp.(stitched_series, 0.0f0, 1.0f0) .* 65535.0f0
    stitched_16 = round.(UInt16, stitched_16)
    save(output_file, stitched_16)
    println("Saved stitched binned series (UInt16) to: $output_file")
    # --- End of your original main function logic ---
end

# Version of main with no arguments - triggers file dialog
function main()
    println("No companion file path provided.")
    println("Opening file selection dialog to choose an .ome.companion file...")

    # Define the filter list for the file dialog
    # Format: "Description1,ext1;ext2;Description2,ext3" or simpler "ext1;ext2"
    filter_list = "companion.ome" # Simple filter for files ending in .companion.ome
    # Alternatively, a more descriptive filter (might vary slightly by OS):
    # filter_list = "OME Companion Files,companion.ome;All Files,*"

    # Use pick_file to open the dialog. Start in the current working directory (pwd()).
    selected_file = pick_file(pwd(), filterlist=filter_list)

    if isempty(selected_file)
        println("File selection cancelled. Exiting.")
        return # Exit gracefully if the user cancels
    else
        println("File selected: ", selected_file)
        # Call the other main function with the selected path
        main(selected_file)
    end
end

#base_path = raw"C:\Users\uComp\Documents\LinneaData\20250226\Sample5_1\\"
#companion_file = base_path * "20250226_ll_sample5_60x_10msexp_1-2bint_gfp-bypass_60sint1.companion.ome"

#base_path = raw"C:\Users\uComp\Downloads\sample\\"
#companion_file = base_path * "20231113_1101_129_e14_kidney_w1_c_dapi_dba_ck8_sox9_10x2.companion.ome"
#main(companion_file)


main()