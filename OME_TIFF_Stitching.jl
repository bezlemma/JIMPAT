using FileIO, EzXML, FFTW, Base.Threads, LinearAlgebra, TiffImages
using Statistics, SparseArrays, OrderedCollections
using NativeFileDialog, ProgressMeter
using GLMakie

# Assuming these files contain the necessary function definitions
include("Stitching_Load_Functions.jl")
include("Stitching_Fuse_Functions.jl")
include("Stitching_Align_Functions.jl")
include("SaveTiffPlz.jl")

# --- Helper Function for Parsing Range Input (Modified for GUI) ---

"""
    parse_range_input_gui(input_str, max_val, dim_name)

Parses user input string for a range (e.g., "1:10", "5", "").
Returns a valid UnitRange or throws an informative error message string.
"""
function parse_range_input_gui(input_str::AbstractString, max_val::Int, dim_name::String)
    input_str = strip(input_str)
    default_range = 1:max_val

    if isempty(input_str)
        return default_range # Return default if empty
    end

    try
        # Try parsing single integer
        if occursin(r"^\d+$", input_str)
            val = parse(Int, input_str)
            if 1 <= val <= max_val
                return val:val
            else
                # Return error message instead of throwing
                return "Single value $val out of valid range [1, $max_val] for $dim_name."
            end
        end

        # Try parsing start:end format
        if occursin(r"^\d+:\d+$", input_str)
            parts = split(input_str, ':')
            start_val = parse(Int, parts[1])
            end_val = parse(Int, parts[2])

            if start_val > end_val
                return "Start value $start_val cannot be greater than end value $end_val for $dim_name."
            end
            if 1 <= start_val <= max_val && 1 <= end_val <= max_val
                return start_val:end_val
            else
                return "Range $start_val:$end_val contains values out of valid range [1, $max_val] for $dim_name."
            end
        end

        # If neither format matches
        return "Invalid format for $dim_name range: '$input_str'. Use 'start:end', a single number, or leave blank for default (1:$max_val)."

    catch e # Catch parsing errors (e.g., integer overflow)
        return "Error parsing $dim_name range input: $(sprint(showerror, e))"
    end
end


# --- GUI Function for Range Input ---

"""
    prompt_for_ranges_gui(max_t, max_z, max_c)

Opens a GUI window using LibUI.jl to prompt the user for T, Z, C ranges.
Returns a tuple `(t_range, z_range, c_range)` if OK is clicked and parsing succeeds,
or `nothing` if Cancel is clicked or the window is closed.
"""
function prompt_for_ranges_gui(max_t::Int, max_z::Int, max_c::Int)
    # Channel to signal completion and pass back results from the GUI callback
    result_channel = Channel{Union{Tuple{UnitRange,UnitRange,UnitRange},Nothing}}(1)

    # Initialize LibUI event loop on the main thread if not already running
    # LibUI.jl handles this internally now, usually don't need manual init/quit
    # ui_task = @async LibUI.main() # Run event loop in background task if needed

    win = LibUI.Window("Select Processing Ranges", 300, 150, false)
    LibUI.set_margined(win, true)
    LibUI.on_closing(win) do w
        println("GUI Window closed by user.")
        try
            # Check if channel is still open before trying to put!
            # Check isempty to ensure we dont overwrite a value already put by OK/Cancel
            if isopen(result_channel) && isempty(result_channel)
                put!(result_channel, nothing) # Signal cancellation
            end
        catch e
            @warn "Error putting nothing to channel on window close: $e"
        end
        LibUI.destroy(win) # Ensure window resources are freed
        # If LibUI.main() was started manually, might need LibUI.quit() here or after take!
        return true # Allow closing
    end

    vbox = LibUI.VerticalBox()
    LibUI.set_padded(vbox, true)

    # --- T Range Input ---
    hbox_t = LibUI.HorizontalBox()
    LibUI.set_padded(hbox_t, true)
    label_t = LibUI.Label("Time (T) Range [1:$max_t]:")
    entry_t = LibUI.Entry()
    LibUI.set_text(entry_t, "1:$max_t") # Default text
    LibUI.append(hbox_t, label_t, false)
    LibUI.append(hbox_t, entry_t, true)
    LibUI.append(vbox, hbox_t, false)

    # --- Z Range Input ---
    hbox_z = LibUI.HorizontalBox()
    LibUI.set_padded(hbox_z, true)
    label_z = LibUI.Label("Z Range [1:$max_z]:   ") # Pad for alignment
    entry_z = LibUI.Entry()
    LibUI.set_text(entry_z, "1:$max_z") # Default text
    LibUI.append(hbox_z, label_z, false)
    LibUI.append(hbox_z, entry_z, true)
    LibUI.append(vbox, hbox_z, false)

    # --- C Range Input ---
    hbox_c = LibUI.HorizontalBox()
    LibUI.set_padded(hbox_c, true)
    label_c = LibUI.Label("Channel (C) Range [1:$max_c]:")
    entry_c = LibUI.Entry()
    LibUI.set_text(entry_c, "1:$max_c") # Default text
    LibUI.append(hbox_c, label_c, false)
    LibUI.append(hbox_c, entry_c, true)
    LibUI.append(vbox, hbox_c, false)

    # --- Buttons ---
    hbox_buttons = LibUI.HorizontalBox()
    LibUI.set_padded(hbox_buttons, true)

    ok_button = LibUI.Button("OK")
    cancel_button = LibUI.Button("Cancel")

    LibUI.append(hbox_buttons, ok_button, true)
    LibUI.append(hbox_buttons, cancel_button, true)
    LibUI.append(vbox, hbox_buttons, false)

    LibUI.set_child(win, vbox)

    # --- Button Callbacks ---
    LibUI.on_clicked(ok_button) do btn
        t_text = LibUI.text(entry_t)
        z_text = LibUI.text(entry_z)
        c_text = LibUI.text(entry_c)

        t_result = parse_range_input_gui(t_text, max_t, "Time (T)")
        z_result = parse_range_input_gui(z_text, max_z, "Z")
        c_result = parse_range_input_gui(c_text, max_c, "Channel (C)")

        # Check if any parsing failed (returned an error string)
        errors = filter(x -> x isa String, [t_result, z_result, c_result])

        if !isempty(errors)
            error_message = join(errors, "\n")
            println("Input errors:\n$error_message")
            # Display error message in a GUI popup
            LibUI.msgbox_error(win, "Input Error", error_message)
        else
            # All results are valid ranges
            println("Ranges accepted: T=$(t_result), Z=$(z_result), C=$(c_result)")
            try
                if isopen(result_channel) && isempty(result_channel)
                    put!(result_channel, (t_result, z_result, c_result))
                end
            catch e
                @warn "Error putting result tuple to channel on OK: $e"
            end
            LibUI.destroy(win) # Close the window
            # If LibUI.main() was started manually, might need LibUI.quit()
        end
    end

    LibUI.on_clicked(cancel_button) do btn
        println("Cancel button clicked.")
        try
            if isopen(result_channel) && isempty(result_channel)
                put!(result_channel, nothing) # Signal cancellation
            end
        catch e
            @warn "Error putting nothing to channel on Cancel: $e"
        end
        LibUI.destroy(win) # Close the window
        # If LibUI.main() was started manually, might need LibUI.quit()
    end

    LibUI.show(win)

    # Block and wait for the result from the channel
    # Use timedwait or similar if you want to prevent indefinite blocking
    println("Waiting for GUI input...")
    final_result = take!(result_channel)
    close(result_channel)
    println("GUI interaction complete.")

    # If we started LibUI.main() manually and need to stop it:
    # LibUI.quit()
    # wait(ui_task) # Wait for the LibUI task to finish if started manually

    return final_result # Return tuple of ranges or nothing
end


# --- Main Execution Logic (Core Function with Ranges - UNCHANGED) ---
# This function remains the same as in the previous version
function main(companion_file::String, time_range::AbstractRange, z_range::AbstractRange, c_range::AbstractRange)
    println("--- Starting Stitching Process ---")
    println("Using Companion File: ", companion_file)
    println("Specified Time Range (T): ", time_range)
    println("Specified Z Range (Z):    ", z_range)
    println("Specified Channel Range (C):", c_range)

    base_path = dirname(companion_file)

    println("\nParsing OME companion file...")
    # ... (rest of the parsing logic remains the same)
    image_data = parse_ome_companion(companion_file)
    if isempty(image_data)
        println("Error: Parsing failed or returned no valid tile data.")
        return
    end
    println("Parsed data for $(length(image_data)) valid tiles.")
    if !haskey(image_data[1], "metadata_size_c") || !haskey(image_data[1], "metadata_size_z") || !haskey(image_data[1], "metadata_size_t")
        println("Error: Cannot determine metadata dimensions from the first tile's data.")
        return
    end
    metadata_size_t = image_data[1]["metadata_size_t"]
    metadata_size_z = image_data[1]["metadata_size_z"]
    metadata_size_c = image_data[1]["metadata_size_c"]
    println("Dimensions from Metadata: T=$metadata_size_t, Z=$metadata_size_z, C=$metadata_size_c")
    if metadata_size_t <= 0 || metadata_size_z <= 0 || metadata_size_c <= 0
        println("Error: Metadata dimensions reported as zero or negative.")
        return
    end

    # --- Validate User-Provided Ranges ---
    if maximum(time_range) > metadata_size_t || minimum(time_range) < 1
        println("Error: Provided Time range ($(time_range)) is outside the valid metadata range (1:$(metadata_size_t)).")
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
    println("Provided ranges are valid.")

    # --- Parameters ---
    align_strategy = "From Position"
    bin_factor = 1
    Fusion_Strategy = "None"
    feather_blend_size = 30.0

    offsets = tile_using_positions(image_data, bin_factor)
    println("Calculated offsets for $(length(offsets)) tiles.")

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
    println("--- Starting Stitching Process (Default Ranges) ---")
    println("Using Companion File: ", companion_file)

    println("\nParsing OME companion file to determine dimensions...")
    # ... (Parsing logic to get metadata sizes remains the same)
    image_data = parse_ome_companion(companion_file)
    if isempty(image_data) # ... (error checks remain the same)
        println("Error: Parsing failed or returned no valid tile data.")
        return
    end
    if !haskey(image_data[1], "metadata_size_c") # ... (error checks remain the same)
        println("Error: Cannot determine metadata dimensions from the first tile's data.")
        return
    end
    metadata_size_t = image_data[1]["metadata_size_t"]
    metadata_size_z = image_data[1]["metadata_size_z"]
    metadata_size_c = image_data[1]["metadata_size_c"]
    println("Dimensions from Metadata: T=$metadata_size_t, Z=$metadata_size_z, C=$metadata_size_c")
    if metadata_size_t <= 0 || metadata_size_z <= 0 || metadata_size_c <= 0 # ... (error checks remain the same)
        println("Error: Metadata dimensions reported as zero or negative.")
        return
    end

    default_time_range = 1:metadata_size_t
    default_z_range = 1:metadata_size_z
    default_c_range = 1:metadata_size_c
    println("Using default full ranges derived from metadata.")

    main(companion_file, default_time_range, default_z_range, default_c_range)
end

"""
    prompt_for_ranges_makie(max_t, max_z, max_c)
Opens a GUI window using GLMakie to prompt the user for T, Z, C ranges.
Returns a tuple `(t_range, z_range, c_range)` if OK is clicked and parsing succeeds,
or `nothing` if Cancel/Close is clicked.
"""
function prompt_for_ranges_makie(max_t::Int, max_z::Int, max_c::Int)
    result_channel = Channel{Union{Tuple{UnitRange,UnitRange,UnitRange},Nothing}}(1)

    # Set up the Makie figure and layout
    fig = Figure(size=(400, 250), title="Select Processing Ranges")
    gl = fig[1, 1] = GridLayout()

    # Add Labels and Textboxes using Makie widgets
    label_t = Label(gl[1, 1], "Time (T) Range [1:$max_t]:", halign=:right)
    box_t = Textbox(gl[1, 2], placeholder="1:$max_t", width=150) # Use placeholder or set_text
    box_t.stored_string[] = "1:$max_t" # Set initial text

    label_z = Label(gl[2, 1], "Z Range [1:$max_z]:", halign=:right)
    box_z = Textbox(gl[2, 2], placeholder="1:$max_z", width=150)
    box_z.stored_string[] = "1:$max_z"

    label_c = Label(gl[3, 1], "Channel (C) Range [1:$max_c]:", halign=:right)
    box_c = Textbox(gl[3, 2], placeholder="1:$max_c", width=150)
    box_c.stored_string[] = "1:$max_c"

    # Error message label (initially empty)
    error_label = Label(gl[4, 1:2], "", color=:red, tellwidth=false)

    # OK and Cancel Buttons
    button_layout = gl[5, 1:2] = GridLayout()
    ok_button = Button(button_layout[1, 1], label="OK", width=100)
    cancel_button = Button(button_layout[1, 2], label="Cancel", width=100)

    # --- Callbacks ---
    on(ok_button.clicks) do _
        t_text = box_t.stored_string[]
        z_text = box_z.stored_string[]
        c_text = box_c.stored_string[]

        t_result = parse_range_input_gui(t_text, max_t, "Time (T)")
        z_result = parse_range_input_gui(z_text, max_z, "Z")
        c_result = parse_range_input_gui(c_text, max_c, "Channel (C)")

        errors = filter(x -> x isa String, [t_result, z_result, c_result])

        if !isempty(errors)
            error_message = join(errors, "\n")
            error_label.text = error_message # Display error in the GUI label
            println("Input errors:\n$error_message") # Also print to console
        else
            # Success!
            error_label.text = "" # Clear any previous error message
            println("Ranges accepted: T=$(t_result), Z=$(z_result), C=$(c_result)")
            try
                # Ensure channel is open and empty before putting result
                if isopen(result_channel) && isempty(result_channel)
                    put!(result_channel, (t_result, z_result, c_result))
                end
                # Close the window - get the screen and close it
                screen = GLMakie.Screen(fig.scene)
                close(screen)
            catch e
                @warn "Error putting result or closing window on OK: $e"
            end
        end
    end

    on(cancel_button.clicks) do _
        println("Cancel button clicked.")
        try
            if isopen(result_channel) && isempty(result_channel)
                put!(result_channel, nothing)
            end
            screen = GLMakie.Screen(fig.scene)
            close(screen)
        catch e
            @warn "Error putting nothing or closing window on Cancel: $e"
        end
    end

    # Handle window closing via 'X' button
    # Makie handles window closing events internally; need to ensure the channel gets 'nothing'
    # One way is to check the channel state *after* the display call returns
    # Or rely on a timeout for take! below

    # Display the figure - this opens the window
    # The `display` function might block or might return immediately depending on context.
    # Running it within a script often keeps the window open until script ends or window closed.
    println("Displaying Makie window...")
    screen = display(fig) # Or display(GLMakie.Screen(), fig)

    # Block and wait for the result from the channel
    # Add a timeout (e.g., Inf or a large number) to wait indefinitely
    # Or handle the case where the window is closed manually via 'X'
    final_result = nothing
    try
        println("Waiting for GUI input...")
        # Wait indefinitely until something is put on the channel
        final_result = take!(result_channel)
        println("GUI interaction complete (OK/Cancel).")
    catch e
        println("Error taking from channel (potentially closed early?): $e")
        # Ensure window is closed if take! fails or is interrupted
        if isopen(screen)
            close(screen)
        end
        final_result = nothing # Assume cancellation if error
    end

    # Check if window was closed manually (if take! didn't get a value)
    if isempty(result_channel) && isnothing(final_result)
        println("Window closed manually without OK/Cancel.")
        final_result = nothing # Ensure it's nothing
    end

    close(result_channel) # Close the channel

    # Ensure the screen is definitely closed if it somehow wasn't
    # Although callbacks should handle this now
    if isopen(screen)
        println("Closing lingering screen...")
        try
            close(screen)
        catch
        end
    end

    return final_result
end

# --- main() function (no arguments) - MODIFIED to use Makie GUI ---
function main()
    println("--- Starting Stitching Process (Interactive GUI Mode - Makie) ---")
    # ... (File selection using NativeFileDialog remains the same) ...
    println("Opening file selection dialog...")
    filter_list = "companion.ome"
    selected_file = ""
    try
        selected_file = pick_file(pwd(); filterlist=filter_list)
        if isempty(selected_file)
            println("File selection cancelled. Exiting.")
            return
        end
        println("File selected: ", selected_file)
    catch e
        println("\nError during file selection dialog:")
        showerror(stdout, e)
        println("\nExiting.")
        return
    end

    # ... (Parsing metadata remains the same) ...
    println("\nParsing selected file to determine dimensions...")
    local image_data
    try # ... (parsing logic) ...
        image_data = parse_ome_companion(selected_file)
        if isempty(image_data) #... (error check) ...
            println("Error: Parsing failed or returned no valid tile data from selected file.")
            return
        end
        if !haskey(image_data[1], "metadata_size_c") # ... (error check) ...
            println("Error: Cannot determine metadata dimensions from the first tile's data.")
            return
        end
    catch e # ... (error check) ...
        println("\nError parsing the selected companion file:")
        showerror(stdout, e)
        println("\nExiting.")
        return
    end
    metadata_size_t = image_data[1]["metadata_size_t"]
    metadata_size_z = image_data[1]["metadata_size_z"]
    metadata_size_c = image_data[1]["metadata_size_c"]
    println("Dimensions from Metadata: T=$metadata_size_t, Z=$metadata_size_z, C=$metadata_size_c")
    if metadata_size_t <= 0 || metadata_size_z <= 0 || metadata_size_c <= 0 # ... (error check) ...
        println("Error: Metadata dimensions reported as zero or negative.")
        return
    end

    # --- Prompt for ranges using the Makie GUI ---
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



