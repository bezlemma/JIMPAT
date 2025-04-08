module Stitching

using FileIO, Printf, FFTW, Base.Threads, LinearAlgebra, TiffImages
using Statistics, SparseArrays, OrderedCollections
using NativeFileDialog, ProgressMeter
using EzXML, DataFrames
using ArgParse

export stitch # Export the primary function

include("Stitching_Load_Functions.jl")
include("Stitching_Fuse_Functions.jl")
include("Stitching_Align_Functions.jl")
#include("Stitching_GUI.jl")
include("SaveTiffPlz.jl")

# --- (Optional but recommended) Keep internal core logic separate ---
# This function does the actual work once all parameters are known
function _execute_stitching(companion_file::String, t_range::UnitRange{Int}, z_range::UnitRange{Int}, c_range::UnitRange{Int},
    save_location::String, bin_factor::Int, align_strategy::String, fusion_strategy::String, feather_blend_size::Float64)
    #Start by giving back user their information
    println("--- Executing Stitching ---")
    println("  Input File:        ", companion_file); println("  Save Location:     ", save_location);
    println("  Time Range (T):    ", t_range); println("  Z Range:           ", z_range); println("  Channel Range (C): ", c_range);
    println("  Bin Factor:        ", bin_factor); println("  Align Strategy:    ", align_strategy); println("  Fusion Strategy:   ", fusion_strategy);
    if fusion_strategy == "Feather"; println("  Feather Blend Size: ", feather_blend_size); end;
    println("--------------------------")

    # Algin and Fuse
    base_path = dirname(companion_file)
    image_data = parse_ome(companion_file)
    println("\nStarting alignment...")
    offsets = tile_using_positions(image_data, bin_factor) # Assuming only this strategy for now
    println("\nStarting fusion...")
    stitched_series = fuse_main(image_data, offsets, fusion_strategy, base_path,
        c_range, z_range, t_range, bin_factor, feather_blend_size)

    # Save
    mkpath(save_location)
    binned_str = bin_factor > 1 ? "_binned$(bin_factor)x" : ""
    fusion_str = fusion_strategy != "None" ? "_$(lowercase(fusion_strategy))" : ""
    output_filename = replace(basename(companion_file), r"(\.companion)?\.ome$"i => "") * "_stitched$(binned_str)$(fusion_str).tif"
    output_file = joinpath(save_location, output_filename)
    println("Saving to: $output_file")
    stitched_16 = round.(UInt16, clamp.(Float32.(stitched_series), 0.0f0, 1.0f0) .* 65535.0f0)
    SaveTiffPlz(output_file, stitched_16)
    println("Saving complete.")
end

# --- Primary User-Facing Function (Handles Scripting & Internal Calls) ---
"""
    stitch(; loadfile=nothing, t_range=nothing, ..., kwargs...)
Runs the stitching process. Handles file picking via GUI if `loadfile`
is not provided, sets defaults, resolves parameters, and executes stitching.
"""
function stitch(; loadfile::Union{String, Nothing}=nothing,
        t_range::Union{AbstractString, AbstractRange, Integer, Nothing} = nothing,
        z_range::Union{AbstractString, AbstractRange, Integer, Nothing} = nothing,
        c_range::Union{AbstractString, AbstractRange, Integer, Nothing} = nothing,
        save_location::Union{String, Nothing} = nothing,
        bin_factor::Int = 4,
        align_strategy::String = "From Position",
        fusion_strategy::String = "None",
        feather_blend_size::Float64 = 30.0
    )

    # --- Handle File Input (GUI if needed) ---
    actual_loadfile = loadfile # Start with provided value
    if isnothing(actual_loadfile) || (isa(actual_loadfile, AbstractString) && isempty(actual_loadfile))
        println("Input file not specified, opening GUI file picker...")
            filter_list = "ome"
            actual_loadfile = pick_file(pwd(); filterlist=filter_list)
            if isempty(actual_loadfile)
                error("No file selected from GUI. Aborting.")
            end
    end

    # --- Validate Final File Path ---
    if !isfile(actual_loadfile)
        error("Input file not found or invalid: $actual_loadfile")
    end
    println("Using Input File: ", actual_loadfile)

    # --- Parse Metadata ---
    println("Parsing metadata...")
    local meta_data
    try
        meta_data = parse_ome_metadata(actual_loadfile) # Assumes parse_ome_metadata is defined
    catch e
        println("ERROR: Failed to parse metadata from '$actual_loadfile'.")
        rethrow(e)
    end
    metadata_size_t = meta_data["SizeT"]
    metadata_size_z = meta_data["SizeZ"]
    metadata_size_c = meta_data["SizeC"]
    println("Metadata dimensions: T=$(metadata_size_t), Z=$(metadata_size_z), C=$(metadata_size_c)")

    function _parse_range_str(range_str, max_dim)
        range_str_l = lowercase(strip(range_str))
        if range_str_l == "all" return 1:max_dim end
        if occursin(":", range_str_l)
            parts = split(range_str_l, ':'); start, stop = parse.(Int, parts)
        else
            start = stop = parse(Int, range_str_l)
        end
        if !(1 <= start <= stop <= max_dim) error("Range $start:$stop out of bounds (1:$max_dim)") end
        return start:stop
    end

    function _resolve_range(val, max_dim, dim_name)
        if isnothing(val) return 1:max_dim end # Default = all
        try
            if isa(val, AbstractString) return _parse_range_str(val, max_dim) end
            if isa(val, Integer) val_int=Int(val); return val_int:val_int end
            if isa(val, AbstractRange) min_v,max_v=minimum(val),maximum(val); return UnitRange{Int}(min_v:max_v) end
        catch e
             error("Invalid format for $dim_name range '$val'. Use 'all', 'N', 'N:M', Int, or Range. Error: $e")
        end
        error("Unhandled type for $dim_name range: $(typeof(val))") # Failsafe
    end

    local actual_t_range, actual_z_range, actual_c_range
    try
        actual_t_range = _resolve_range(t_range, metadata_size_t, "T")
        actual_z_range = _resolve_range(z_range, metadata_size_z, "Z")
        actual_c_range = _resolve_range(c_range, metadata_size_c, "C")
    catch e; rethrow(e); end

    # Resolve Save Location (default = input file dir)
    actual_save_location = isnothing(save_location) ? dirname(actual_loadfile) : save_location
    if !isa(actual_save_location, AbstractString) || isempty(actual_save_location)
         error("Invalid save_location derived.")
    end

    # --- Call Core Logic ---
    _execute_stitching(
        actual_loadfile,
        actual_t_range, actual_z_range, actual_c_range,
        actual_save_location, bin_factor, align_strategy,
        fusion_strategy, feather_blend_size
    )
end


# --- Command-Line Argument Parsing Setup ---
function setup_argparse()
    s = ArgParseSettings(description="Stitch tiled images from OME.", prog="Stitching.jl")
    @add_arg_table! s begin
        "--loadfile", "-l"
            help = "Path to the .companion.ome file. If omitted, a GUI picker appears."
            arg_type = String; default = nothing;
        "--trange", "-t"
            help = "Time range ('all', 'N', 'N:M'). Defaults to 'all'."
            arg_type = String; default = "all";
        "--zrange", "-z"
            help = "Z range ('all', 'N', 'N:M'). Defaults to 'all'."
            arg_type = String; default = "all";
        "--crange", "-c"
            help = "Channel range ('all', 'N', 'N:M'). Defaults to 'all'."
            arg_type = String; default = "all";
        "--savelocation", "-s"
            help = "Save directory. Defaults to input file directory."
            arg_type = String; default = nothing; 
        "--binsize", "-b"
            help = "Binning factor."
            arg_type = Int; default = 4;
        "--align_strategy", "-a"
            help = "Alignment strategy ('From Position')."
            arg_type = String; default = "From Position";
        "--fusion_strategy", "-f"
            help = "Fusion strategy ('None', 'Feather')."
            arg_type = String; default = "None";
         "--feather_blend_size"
             help = "Blend size for Feather fusion."
             arg_type = Float64; default = 30.0;
    end
    return s
end

# --- Main Entry Point for Command-Line Execution ---
function julia_main()::Cint
    arg_settings = setup_argparse()
    local parsed_args
    try
        # Use parse_args(ARGS, settings; as_symbols=true) for direct kwarg passing
        # Requires Julia 1.0+ for as_symbols
        parsed_args = parse_args(ARGS, arg_settings; as_symbols = true)
    catch e
        # ArgParse usually handles errors/help, but catch others
        println("Error parsing arguments: $e")
        return 1
    end

    # Remove the :%COMMAND% entry if present (older ArgParse versions)
    pop!(parsed_args, Symbol("%COMMAND%"), nothing)
    try
        stitch(; parsed_args...) # Splat the dictionary as keyword arguments
        return 0 # Success
    catch e
        println("\n ERROR during stitching")
        return 1 # Failure
    end
end

# Execute main if script is run directly
if abspath(PROGRAM_FILE) == @__FILE__; julia_main(); end;

end