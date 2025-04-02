
# --- Blend Functions ---

"""
Main fusing function orchestrating the process for C, Z, T ranges.
"""
function fuse_main(
    image_data, offsets::Dict{Any,Tuple{Float64,Float64}}, Fusion_Strategy::String,
    base_path::String, c_range, z_range, time_range, bin_factor::Int,
    metadata_size_c::Int, metadata_size_z::Int, metadata_size_t::Int,
    feather_blend_size::Float64)

    num_c, num_z, num_t = length(c_range), length(z_range), length(time_range)
    println("Bin Factor: $bin_factor")
    println("Fusion Strategy: $Fusion_Strategy")
    if Fusion_Strategy == "Feather"
        println("Feather Size: $feather_blend_size")
    end
    if isempty(image_data) || isempty(offsets)
        println("Error: No image data or offsets.")
        return Array{Float32,5}(undef, 0, 0, 0, 0, 0)
    end

    first_c, first_z, first_t = first(c_range), first(z_range), first(time_range)
    println("Stitching reference plane (C=$first_c, Z=$first_z, T=$first_t) to determine size...")
    local mosaic_ref
    local h0, w0
    stitch_func = Fusion_Strategy == "Feather" ? stitch_plane_feather : stitch_plane
    kwargs = Fusion_Strategy == "Feather" ? (; feather_size=feather_blend_size) : (;)
    mosaic_ref = stitch_func(image_data, offsets, base_path; c_index=first_c, z_index=first_z, time_index=first_t, b_factor=bin_factor, kwargs...)
    if isempty(mosaic_ref) || size(mosaic_ref) == (1, 1) && all(mosaic_ref .== 0.0f0)
        println("Error: Ref plane stitching failed.")
        return Array{Float32,5}(undef, 0, 0, 0, 0, 0)
    end
    (h0, w0) = size(mosaic_ref)

    stitched_series = Array{Float32,5}(undef, h0, w0, num_c, num_z, num_t)
    println("Allocated 5D output array: ($h0, $w0, $num_c, $num_z, $num_t)")

    ref_ci = findfirst(==(first_c), c_range)
    ref_zi = findfirst(==(first_z), z_range)
    ref_ti = findfirst(==(first_t), time_range)
    if !isnothing(ref_ci) && !isnothing(ref_zi) && !isnothing(ref_ti)
        stitched_series[:, :, ref_ci, ref_zi, ref_ti] = mosaic_ref
    else
        mosaic_ref = nothing
    end

    @showprogress @threads for idx in CartesianIndices((num_c, num_z, num_t))
        ci, zi, ti = Tuple(idx) # Get 1-based indices for output array

        # Map back to actual C/Z/T values from ranges
        c = c_range[ci]
        z = z_range[zi]
        t = time_range[ti]

        # Skip if this is the reference plane and it was already placed
        if !isnothing(mosaic_ref) && c == first_c && z == first_z && t == first_t
            continue
        end

        try
            local stitched_plane::Matrix{Float32}
            stitch_func = Fusion_Strategy == "Feather" ? stitch_plane_feather : stitch_plane
            kwargs = Fusion_Strategy == "Feather" ? (; feather_size=feather_blend_size) : (;)
            stitched_plane = stitch_func(image_data, offsets, base_path; c_index=c, z_index=z, time_index=t, b_factor=bin_factor, kwargs...)

            if size(stitched_plane) == (h0, w0)
                stitched_series[:, :, ci, zi, ti] = stitched_plane
            else
                println("Warn: Plane C=$c,Z=$z,T=$t size $(size(stitched_plane))!=($h0,$w0). Zeros.")
                stitched_series[:, :, ci, zi, ti] .= 0.0f0
            end
        catch e
            println("\nError stitching plane C=$c, Z=$z, T=$t:")
            showerror(stdout, e)
            println()
            println("Filling plane C=$c, Z=$z, T=$t with zeros.")
            # Ensure assignment happens even on error to avoid uninitialized array elements
            # Check bounds before assigning zeros
            if 1 <= ci <= size(stitched_series, 3) && 1 <= zi <= size(stitched_series, 4) && 1 <= ti <= size(stitched_series, 5)
                try
                    stitched_series[:, :, ci, zi, ti] .= 0.0f0
                catch
                end # Use try-catch for assignment robustness in threads
            end
        end
    end # End of @threads loop

    return stitched_series
end

"""
Stitches a single C, Z, T plane using simple overwrite ("None" fusion).
Finds filename based on parsed metadata for the given channel index.
"""
function stitch_plane(
    image_data, offsets::Dict{Any,Tuple{Float64,Float64}}, base_path::String;
    c_index::Int=1, z_index::Int=1, time_index::Int=1, b_factor::Int=1)

    relevant_tiles = [tile for tile in image_data if haskey(offsets, tile) && haskey(tile, "channel_files")]
    if isempty(relevant_tiles)
        return fill(0.0f0, 1, 1)
    end

    tile_images = Dict{Any,Matrix{Float32}}()
    xvals = Float64[]
    yvals = Float64[]

    for tile in relevant_tiles
        channel_files_key = c_index - 1
        if !haskey(tile["channel_files"], channel_files_key)
            continue
        end
        fname = tile["channel_files"][channel_files_key]
        if ismissing(fname)
            continue
        end
        path = get_file_path(base_path, fname)

        try
            tile_im = load_plane_and_bin(path, c_index, z_index, time_index, b_factor=b_factor)
            if isempty(tile_im)
                continue
            end
            tile_images[tile] = tile_im
            (ox, oy) = offsets[tile]
            h_tile, w_tile = size(tile_im)
            if h_tile > 0 && w_tile > 0
                push!(xvals, ox, ox + w_tile - 1)
                push!(yvals, oy, oy + h_tile - 1)
            else
                println("Warn: Tile $(tile["id"]) Z=$z_index, C=$c_index zero dim. Skip.")
            end
        catch e
            println("Skip tile $(tile["id"]) C=$c_index, Z=$z_index due to error. File: $path")
            continue
        end
    end

    if isempty(xvals)
        @warn "NoneFusion: No valid tiles for C=$c_index, Z=$z_index, T=$time_index."
        return fill(0.0f0, 1, 1)
    end
    min_x = floor(Int, minimum(xvals))
    max_x = ceil(Int, maximum(xvals))
    min_y = floor(Int, minimum(yvals))
    max_y = ceil(Int, maximum(yvals))
    width = max_x - min_x + 1
    height = max_y - min_y + 1
    if width <= 0 || height <= 0
        @warn "NoneFusion: Invalid canvas dims ($height,$width) C=$c_index, Z=$z_index, T=$time_index."
        return fill(0.0f0, 1, 1)
    end

    canvas = fill(0.0f0, height, width)
    for tile in keys(tile_images)
        im = tile_images[tile]
        (ox, oy) = offsets[tile]
        h_tile, w_tile = size(im)
        oy_c_start = round(Int, oy - min_y) + 1
        ox_c_start = round(Int, ox - min_x) + 1
        oy_c_end = oy_c_start + h_tile - 1
        ox_c_end = ox_c_start + w_tile - 1
        y_start = max(1, oy_c_start)
        y_end = min(height, oy_c_end)
        x_start = max(1, ox_c_start)
        x_end = min(width, ox_c_end)
        im_y_start = y_start - oy_c_start + 1
        im_y_end = y_end - oy_c_start + 1
        im_x_start = x_start - ox_c_start + 1
        im_x_end = x_end - ox_c_start + 1
        if y_start <= y_end && x_start <= x_end && im_y_start <= im_y_end && im_x_start <= im_x_end
            cv = view(canvas, y_start:y_end, x_start:x_end)
            iv = view(im, im_y_start:im_y_end, im_x_start:im_x_end)
            if size(cv) == size(iv)
                cv .= iv
            else
                println("Warn: Size mismatch $(tile["id"]). Skip.")
            end
        end
    end
    return canvas
end

"""
Stitches a single C, Z, T plane using Feather blending.
Finds filename based on parsed metadata for the given channel index.
"""
function stitch_plane_feather(
    image_data, offsets::Dict{Any,Tuple{Float64,Float64}}, base_path::String;
    c_index::Int=1, z_index::Int=1, time_index::Int=1, b_factor::Int=1, feather_size::Float64=30.0)

    relevant_tiles = [tile for tile in image_data if haskey(offsets, tile) && haskey(tile, "channel_files")]
    if isempty(relevant_tiles)
        return fill(0.0f0, 1, 1)
    end

    tile_images = Dict{Any,Matrix{Float32}}()
    xvals = Float64[]
    yvals = Float64[]

    for tile in relevant_tiles
        channel_files_key = c_index - 1
        if !haskey(tile["channel_files"], channel_files_key)
            continue
        end
        fname = tile["channel_files"][channel_files_key]
        if ismissing(fname)
            continue
        end
        path = get_file_path(base_path, fname)
        try
            tile_im = load_plane_and_bin(path, c_index, z_index, time_index, b_factor=b_factor)
            if isempty(tile_im)
                continue
            end
            tile_images[tile] = tile_im
            (ox, oy) = offsets[tile]
            h_tile, w_tile = size(tile_im)
            if h_tile > 0 && w_tile > 0
                push!(xvals, ox, ox + w_tile - 1)
                push!(yvals, oy, oy + h_tile - 1)
            else
                println("Warn: Feather Tile $(tile["id"]) Z=$z_index, C=$c_index zero dim. Skip.")
            end
        catch e
            println("Skip tile $(tile["id"]) C=$c_index, Z=$z_index in Feather. File: $path")
            continue
        end
    end

    if isempty(xvals)
        @warn "Feather: No valid tiles loaded C=$c_index, Z=$z_index, T=$time_index."
        return fill(0.0f0, 1, 1)
    end
    min_x = floor(Int, minimum(xvals))
    max_x = ceil(Int, maximum(xvals))
    min_y = floor(Int, minimum(yvals))
    max_y = ceil(Int, maximum(yvals))
    width = max_x - min_x + 1
    height = max_y - min_y + 1
    if width <= 0 || height <= 0
        @warn "Feather: Invalid canvas dims ($height,$width) C=$c_index, Z=$z_index, T=$time_index."
        return fill(0.0f0, 1, 1)
    end

    accum_image = fill(0.0, height, width)
    accum_weight = fill(0.0, height, width)

    for tile in keys(tile_images) # Feather Blending
        im = tile_images[tile]
        (ox, oy) = offsets[tile]
        h_tile, w_tile = size(im)
        mask = build_feather_mask(h_tile, w_tile; feather_size=feather_size)
        oy_c_start = round(Int, oy - min_y) + 1
        ox_c_start = round(Int, ox - min_x) + 1
        oy_c_end = oy_c_start + h_tile - 1
        ox_c_end = ox_c_start + w_tile - 1
        y_start = max(1, oy_c_start)
        y_end = min(height, oy_c_end)
        x_start = max(1, ox_c_start)
        x_end = min(width, ox_c_end)
        if y_start <= y_end && x_start <= x_end
            @inbounds for canvas_j in y_start:y_end, canvas_i in x_start:x_end
                im_j = canvas_j - oy_c_start + 1
                im_i = canvas_i - ox_c_start + 1
                wgt = Float64(mask[im_j, im_i])
                if wgt > 1e-6
                    accum_image[canvas_j, canvas_i] += Float64(im[im_j, im_i]) * wgt
                    accum_weight[canvas_j, canvas_i] += wgt
                end
            end
        end
    end

    canvas = fill(0.0f0, height, width) # Normalize
    @inbounds for i in eachindex(canvas)
        wgt = accum_weight[i]
        if wgt > 1e-6
            canvas[i] = Float32(accum_image[i] / wgt)
        end
    end
    return canvas
end

function build_feather_mask(h::Int, w::Int; feather_size::Float64=30.0)::Matrix{Float32}
    feather_size = max(1.0, feather_size)
    mask = Matrix{Float32}(undef, h, w)
    @inbounds for j in 1:h, i in 1:w
        med = min(Float32(j - 1), Float32(h - j), Float32(i - 1), Float32(w - i))
        mask[j, i] = clamp(med / feather_size, 0.0f0, 1.0f0)
    end
    return mask
end