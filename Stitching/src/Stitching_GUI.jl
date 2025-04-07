function prompt_for_ranges_makie(max_t::Int, max_z::Int, max_c::Int)
    result_channel = Channel{Union{Tuple{UnitRange,UnitRange,UnitRange},Nothing}}(1)
    fig = Figure(size=(400, 250), title="Select Processing Ranges")
    gl = fig[1, 1] = GridLayout()
    label_t = Label(gl[1, 1], "Time (T) Range [1:$max_t]:", halign=:right)
    box_t = Textbox(gl[1, 2], placeholder="1:$max_t", width=150) # Use placeholder or set_text
    box_t.stored_string[] = "1:$max_t" # Set initial text
    label_z = Label(gl[2, 1], "Z Range [1:$max_z]:", halign=:right)
    box_z = Textbox(gl[2, 2], placeholder="1:$max_z", width=150)
    box_z.stored_string[] = "1:$max_z"
    label_c = Label(gl[3, 1], "Channel (C) Range [1:$max_c]:", halign=:right)
    box_c = Textbox(gl[3, 2], placeholder="1:$max_c", width=150)
    box_c.stored_string[] = "1:$max_c"
    error_label = Label(gl[4, 1:2], "", color=:red, tellwidth=false)
    button_layout = gl[5, 1:2] = GridLayout()
    ok_button = Button(button_layout[1, 1], label="OK", width=100)
    cancel_button = Button(button_layout[1, 2], label="Cancel", width=100)
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
    println("Displaying Makie window...")
    screen = display(fig) # Or display(GLMakie.Screen(), fig)

    final_result = nothing
    try
        println("Waiting for GUI input...")
        final_result = take!(result_channel)
        println("GUI interaction complete (OK/Cancel).")
    catch e
        println("Error taking from channel (potentially closed early?): $e")
        if isopen(screen)
            close(screen)
        end
        final_result = nothing # Assume cancellation if error
    end
    if isempty(result_channel) && isnothing(final_result)
        println("Window closed manually without OK/Cancel.")
        final_result = nothing # Ensure it's nothing
    end
    close(result_channel) # Close the channel
    if isopen(screen)
        println("Closing lingering screen...")
        try
            close(screen)
        catch
        end
    end
    return final_result
end