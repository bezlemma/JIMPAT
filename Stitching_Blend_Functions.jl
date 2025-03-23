using ProgressMeter

function stitch_main(image_data, offsets, Fusion_Strategy, base_path,
    time_range, z_range, c_range)

#TODO: Provide a simple average blend
#TODO: Provide an adfanced carving blend

    if Fusion_Strategy == "Feather"
        mosaic_ref = stitch_timepoint_feather(image_data, offsets, 1, base_path; feather_size=30.0)
        (h0, w0) = size(mosaic_ref)
        stitched_series = Array{Float32,3}(undef, h0, w0, lastindex(time_range))
        stitched_series[:, :, 1] = mosaic_ref
    
        @showprogress @threads for t in time_range
            try
                stitched_series[:, :, t] = stitch_timepoint_feather(image_data, offsets, t, base_path; feather_size=30.0)
            catch
                println("$t failed")
            end
        end
    elseif Fusion_Strategy == "None"
            mosaic_ref = stitch_plane(image_data, offsets, base_path; time_index=0, z_index=0, c_index=0, b_factor=bin_factor)
            (h0, w0) = size(mosaic_ref)
            stitched_series = Array{Float32,3}(undef, h0, w0, lastindex(time_range))
            stitched_series[:, :, 1] = mosaic_ref
    
            @showprogress for t in time_range, z in z_range, c in c_range
                try
                    println("$t")
                    stitched_series[:, :, t] = stitch_plane(image_data, offsets, base_path; time_index=t, z_index=z, c_index=c, b_factor=bin_factor)
                catch
                    println("$t failed")
                end
            end
    end
    return stitched_series
end


function build_feather_mask(h::Int, w::Int; feather_size::Float64=30.0)::Matrix{Float32}
    mask = Matrix{Float32}(undef, h, w)
    for j in 1:h
        dist_top = j - 1
        dist_bottom = h - j
        for i in 1:w
            dist_left = i - 1
            dist_right = w - i
            min_edge_dist = min(dist_top, dist_bottom, dist_left, dist_right)
            if min_edge_dist >= feather_size
                mask[j, i] = 1.0f0
            else
                mask[j, i] = Float32(min_edge_dist / feather_size)
            end
        end
    end
    return mask
end

function stitch_timepoint_feather(image_data, offsets, time_index::Int, base_path::String; feather_size::Float64=30.0)
    relevant_tiles = [img for img in image_data if haskey(offsets, img)]
    if isempty(relevant_tiles)
        return fill(Float32(0), 1, 1)
    end

    tile_images = Dict{Any,Matrix{Float32}}()
    xvals = Float64[]
    yvals = Float64[]

    # load & bin each tile
    for tile in relevant_tiles
        if haskey(tile["tiff_files"], time_index)
            fname = tile["tiff_files"][time_index]
            path = get_file_path(base_path, fname)
            binned_img = load_and_bin(path)  # returns Float32 array now
            tile_images[tile] = binned_img

            (ox, oy) = offsets[tile]
            w = size(binned_img, 2)
            h = size(binned_img, 1)
            push!(xvals, ox, ox + w - 1)
            push!(yvals, oy, oy + h - 1)
        end
    end

    if isempty(xvals)
        return fill(Float32(0), 1, 1)
    end

    min_x = floor(Int, minimum(xvals))
    max_x = ceil(Int, maximum(xvals))
    min_y = floor(Int, minimum(yvals))
    max_y = ceil(Int, maximum(yvals))
    width = max_x - min_x + 1
    height = max_y - min_y + 1

    # We'll accumulate weighted sums here and then divide by total weight
    accum_image = fill(0.0f0, height, width)
    accum_weight = fill(0.0f0, height, width)

    # Feather-based blending
    for tile in relevant_tiles
        if !haskey(tile_images, tile)
            continue
        end
        im = tile_images[tile]
        (ox, oy) = offsets[tile]

        ox_s = round(Int, ox - min_x)
        oy_s = round(Int, oy - min_y)

        h_tile = size(im, 1);  w_tile = size(im, 2)

        mask = build_feather_mask(h_tile, w_tile; feather_size=feather_size)

        for j in 1:h_tile, i in 1:w_tile
            wgt = mask[j, i]
            if wgt > 0
                accum_image[oy_s + j, ox_s + i] += (im[j, i] * wgt)
                accum_weight[oy_s + j, ox_s + i] += wgt
            end
        end
    end

    # Convert accumulated floats back to Float32
    canvas = fill(Float32(0), height, width)
    for y in 1:height, x in 1:width
        wgt = accum_weight[y, x]
        if wgt > 0
            val = accum_image[y, x] / wgt
            # Remove the clamp to 65535, keep fractional floats:
            canvas[y, x] = Float32(val)
        end
    end

    return canvas
end

# old version
# function stitch_timepoint(image_data, offsets, time_index::Int, base_path::String)
#     relevant_tiles = [img for img in image_data if haskey(offsets, img)]
#     xvals = Float64[]; yvals = Float64[]

#     tile_images = Dict{Any,Matrix{Float32}}()

#     for tile in relevant_tiles
#         fname = tile["tiff_files"][time_index]
#         path = get_file_path(base_path, fname)
#         tile_im = load_and_bin(path)
#         tile_images[tile] = tile_im
#         (ox, oy) = offsets[tile]
#         h_tile, w_tile = size(tile_im)
#         push!(xvals, ox, ox + w_tile - 1)
#         push!(yvals, oy, oy + h_tile - 1)
#     end
#     min_x = floor(Int, minimum(xvals))
#     max_x = ceil(Int, maximum(xvals))
#     min_y = floor(Int, minimum(yvals))
#     max_y = ceil(Int, maximum(yvals))
#     width = max_x - min_x + 1
#     height = max_y - min_y + 1

#     canvas = fill(Float32(0), height, width)

#     for tile in relevant_tiles
#         im = tile_images[tile]
#         (ox, oy) = offsets[tile]
#         ox_i = round(Int, ox) - min_x + 1
#         oy_i = round(Int, oy) - min_y + 1
#         h_tile, w_tile = size(im)
#         canvas[oy_i:oy_i+h_tile-1, ox_i:ox_i+w_tile-1] = im
#     end
#     return canvas
# end

# new version that we are trying to get to work
function stitch_plane(image_data, offsets, base_path::String; time_index::Int=1, z_index::Int=1,
    c_index::Int=1, b_factor::Int=1)
    # 1) Identify tiles that have an offset
    #    (If your logic needs to exclude some tiles, adapt here)
    relevant_tiles = [tile for tile in image_data if haskey(offsets, tile)]
    if isempty(relevant_tiles)
        @warn "No tiles found in offsets that match"
        return fill(Float32(0), 1, 1)
    end

    # 2) Prepare placeholders
    tile_images = Dict{Any,Matrix{Float32}}()
    xvals = Float64[]
    yvals = Float64[]

    # 3) Load and place each tile
    for tile in relevant_tiles
        # (A) Figure out which file to open
        #     If each tile has a single file for all planes, e.g.:
        fname = tile["tiff_files"][time_index]  
        # or if you store one file for each (t,z,c), do something like:
        # fname = tile["tiff_files"][(time_index, z_index, channel_index)]
        # Adjust to match your actual data structure.

        path = get_file_path(base_path, fname)

        # (B) Load the correct plane from the TIF
        #     e.g. if it's a multi-plane TIF, do:
        tile_im = load_and_bin(path; 
            b_factor=b_factor
        )
        # If your load_and_bin only takes `time_index`, you’ll need
        # to add arguments for z_index, channel_index inside that function.

        tile_images[tile] = tile_im

        # (C) Offsets (ox, oy): same for all planes, unless you store them per-plane
        (ox, oy) = offsets[tile][(0,0,0)]

        # (D) Update bounding box
        h_tile, w_tile = size(tile_im)
        push!(xvals, ox, ox + w_tile - 1)
        push!(yvals, oy, oy + h_tile - 1)
    end

    # 4) If no bounding box found, return trivial
    if isempty(xvals)
        return fill(Float32(0), 1, 1)
    end

    # 5) Compute bounding box
    min_x = floor(Int, minimum(xvals))
    max_x = ceil(Int, maximum(xvals))
    min_y = floor(Int, minimum(yvals))
    max_y = ceil(Int, maximum(yvals))
    width  = max_x - min_x + 1
    height = max_y - min_y + 1

    # 6) Allocate the stitched canvas
    canvas = fill(Float32(0), height, width)

    # 7) Overwrite each tile in the canvas
    for tile in relevant_tiles
        im = tile_images[tile]
        (ox, oy) = offsets[tile][(0,0,0)]

        ox_i = round(Int, ox) - min_x + 1
        oy_i = round(Int, oy) - min_y + 1

        h_tile, w_tile = size(im)
        canvas[oy_i:oy_i+h_tile-1, ox_i:ox_i+w_tile-1] = im
    end

    return canvas
end

