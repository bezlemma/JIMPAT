
"""
Main fusing function orchestrating the process for C, Z, T ranges.
"""
function fuse_main(image_data, offsets, Fusion_Strategy::String, base_path::String,
    c_range, z_range, t_range, bin_factor, feather_size)

    num_c, num_z, num_t = length(c_range), length(z_range), length(t_range)

    #Stitch together things once, to figure out how large the entire image will be.
    local mosaic_ref, h0, w0
    if Fusion_Strategy == "None"
        mosaic_ref = fuse_plane(image_data, offsets, base_path, c_range[1], z_range[1], t_range[1], bin_factor)
    elseif Fusion_Strategy == "Feather"
        mosaic_ref = stitch_feather(image_data, offsets, base_path, c_range[1], z_range[1], t_range[1], bin_factor, feather_size)
    end
    (h0, w0) = size(mosaic_ref)

    stitched_series = Array{Float32,5}(undef, h0, w0, num_c, num_z, num_t)
    println("Allocated 5D output array: ($h0, $w0, $num_c, $num_z, $num_t)")

    @showprogress @threads for idx in CartesianIndices((num_c, num_z, num_t))
        ci, zi, ti = Tuple(idx) # Get 1-based indices for output array
        c = c_range[ci]
        z = z_range[zi]
        t = t_range[ti]

        if Fusion_Strategy == "None"
            stitched_series[:, :, ci, zi, ti] = fuse_plane(image_data, offsets, base_path, c, z, t, bin_factor)
        elseif Fusion_Strategy == "Feather"
            stitched_series[:, :, ci, zi, ti] = stitch_feather(image_data, offsets, base_path, c, z, t, bin_factor, feather_size)
        end

    end # End of @threads loop

    return stitched_series
end



function fuse_plane(image_data::DataFrame, offsets, base_path::String, c::Int, z::Int, t::Int, bin_factor)
    filtered = filter(row -> row.c == c && row.z == z && row.t == t, image_data)
    sort!(filtered, :pos_id)

    tile_images = Matrix{Float32}[]
    xvals = Float64[]
    yvals = Float64[]

    for row in eachrow(filtered)
        im_raw = load_image(image_data, base_path, row.pos_id, c, z, t)
        img = Float32.(im_raw)
        if bin_factor == 2
            img = bin2(img)
        elseif bin_factor == 4
            img = bin4(img)
        end
        push!(tile_images, img)
        ox, oy = offsets[row.pos_id]
        h, w = size(img)
        push!(xvals, ox, ox + w - 1)
        push!(yvals, oy, oy + h - 1)
    end

    min_x = floor(Int, minimum(xvals))
    max_x = ceil(Int, maximum(xvals))
    min_y = floor(Int, minimum(yvals))
    max_y = ceil(Int, maximum(yvals))
    width = max_x - min_x + 1
    height = max_y - min_y + 1
    canvas = fill(0.0f0, height, width)

    for (i, im) in enumerate(tile_images)
        row = filtered[i, :]
        ox, oy = offsets[row.pos_id]
        h, w = size(im)
        oy_c_start = round(Int, oy - min_y) + 1
        ox_c_start = round(Int, ox - min_x) + 1
        oy_c_end = oy_c_start + h - 1
        ox_c_end = ox_c_start + w - 1
        y_start = max(1, oy_c_start)
        y_end = min(height, oy_c_end)
        x_start = max(1, ox_c_start)
        x_end = min(width, ox_c_end)
        im_y_start = y_start - oy_c_start + 1
        im_y_end = y_end - oy_c_start + 1
        im_x_start = x_start - ox_c_start + 1
        im_x_end = x_end - ox_c_start + 1
        if y_start ≤ y_end && x_start ≤ x_end &&
           im_y_start ≤ im_y_end && im_x_start ≤ im_x_end
            view(canvas, y_start:y_end, x_start:x_end) .=
                view(im, im_y_start:im_y_end, im_x_start:im_x_end)
        end
    end

    return canvas
end



"""
Stitches a single C, Z, T plane using Feather blending.
Finds filename based on parsed metadata for the given channel index.
"""
function stitch_plane_feather(image_data, offsets, base_path::String, c_index, z_index, time_index, b_factor, feather_size)

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
            tile_im = load_plane(path, c_index, z_index, time_index, b_factor=b_factor)
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