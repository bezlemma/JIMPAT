function tile_using_positions(image_data, bin_factor)
    time_index=0
    valid_tiles = [
        tile for tile in image_data if
        tile["grid_row"] !== missing &&
        tile["grid_col"] !== missing &&
        haskey(tile["tiff_files"], time_index)
    ]

    offsets = Dict{Any,Tuple{Float64,Float64}}()
    for tile in valid_tiles
        px = (tile["position_x"] === missing) ? 0.0 : tile["position_x"]
        py = (tile["position_y"] === missing) ? 0.0 : tile["position_y"]
        # Adjust for bin_factor
        offsets[tile] = (px / bin_factor, py / bin_factor)
    end

    return offsets
end

function align_tiles_reference_timepoint(image_data, base_path;
    time_index=0, overlap_frac=0.1,
    max_allowed_shift=100.0)
    valid_tiles = [img for img in image_data if
                   img["grid_row"] !== missing &&
                   img["grid_col"] !== missing &&
                   haskey(img["tiff_files"], time_index)
    ]
    sort!(valid_tiles, by=x -> (x["grid_row"], x["grid_col"]))

    offsets = Dict{Any,Tuple{Float64,Float64}}()
    tile_images = Dict{Any,Matrix{UInt16}}()

    # load & bin each tile
    for tile in valid_tiles
        fname = tile["tiff_files"][time_index]
        path = get_file_path(base_path, fname)
        tile_images[tile] = load_and_bin(path)
    end

    # naive row/col offset guess
    for tile in valid_tiles
        im = tile_images[tile]
        (h, w) = size(im)
        r = tile["grid_row"]
        c = tile["grid_col"]
        ox = c * (w * (1 - overlap_frac))
        oy = r * (h * (1 - overlap_frac))
        offsets[tile] = (ox, oy)
    end

    println("Refining offsets via local cross-correlation at time_index=$time_index ...")
    # local multi-neighbor refinement:
    refined_offsets = refine_offsets(valid_tiles, tile_images, offsets;
        overlap_frac=overlap_frac,
        max_allowed_shift=max_allowed_shift)

    # Quick row/col montage for sanity check
    all_rows = unique(t["grid_row"] for t in valid_tiles)
    sort!(all_rows)
    row_blocks = []
    for r in all_rows
        row_tiles = filter(t -> t["grid_row"] == r, valid_tiles)
        sort!(row_tiles, by=x -> x["grid_col"])
        row_imgs = [tile_images[t] for t in row_tiles]
        push!(row_blocks, hcat(row_imgs...))
    end
    rowcol_image = vcat(row_blocks...)

    return refined_offsets, rowcol_image
end


function overlap_subimages(imA::Matrix{UInt16}, imB::Matrix{UInt16};
    orientation::Symbol, frac::Float64=0.5)
    hA, wA = size(imA)
    hB, wB = size(imB)

    if orientation == :horizontal
        overlap_widthA = round(Int, frac * wA)
        overlap_widthB = round(Int, frac * wB)
        subA = imA[:, wA-overlap_widthA+1:wA]
        subB = imB[:, 1:overlap_widthB]
        offsetA = (0, wA - overlap_widthA)
        offsetB = (0, 0)
        return subA, subB, offsetA, offsetB

    elseif orientation == :vertical
        overlap_heightA = round(Int, frac * hA)
        overlap_heightB = round(Int, frac * hB)
        subA = imA[hA-overlap_heightA+1:hA, :]
        subB = imB[1:overlap_heightB, :]
        offsetA = (hA - overlap_heightA, 0)
        offsetB = (0, 0)
        return subA, subB, offsetA, offsetB

    else
        error("Unknown orientation $orientation")
    end
end

# ---------------------- find_pairwise_offset ----------------------
function find_pairwise_offset(imA::Matrix{UInt16}, imB::Matrix{UInt16})
    A = Float32.(imA)
    B = Float32.(imB)

    FA = fft(A)
    FB = fft(B)
    R = FA .* conj(FB)
    R ./= abs.(R) .+ 1.0f0
    c = real.(ifft(R))

    peak_idx = argmax(c)
    peak_y, peak_x = Tuple(CartesianIndices(c)[peak_idx])
    Ny, Nx = size(c)
    shift_x = peak_x - 1
    shift_y = peak_y - 1

    if shift_x > Nx ÷ 2
        shift_x -= Nx
    end
    if shift_y > Ny ÷ 2
        shift_y -= Ny
    end

    return (shift_x, shift_y)
end

function refine_offsets(valid_tiles, tile_images, initial_offsets;
    overlap_frac::Float64,
    max_allowed_shift::Float64)

    refined = copy(initial_offsets)

    function find_tile_by_rowcol(r, c)
        for t in valid_tiles
            if t["grid_row"] == r && t["grid_col"] == c
                return t
            end
        end
        return nothing
    end

    for tile in valid_tiles
        r = tile["grid_row"]
        c = tile["grid_col"]

        (old_ox, old_oy) = refined[tile]
        neighbor_positions = [(r, c - 1), (r, c + 1), (r - 1, c), (r + 1, c)]
        local_shifts_x = Float64[]
        local_shifts_y = Float64[]

        for (nr, nc) in neighbor_positions
            neighbor_tile = find_tile_by_rowcol(nr, nc)
            if neighbor_tile === nothing
                continue
            end

            orientation = (nr == r) ? :horizontal : :vertical

            imA = tile_images[neighbor_tile]
            imB = tile_images[tile]

            (ox_neighbor, oy_neighbor) = refined[neighbor_tile]
            (ox_guess, oy_guess) = refined[tile]

            subA, subB, offsetA_inA, offsetB_inB = overlap_subimages(
                imA, imB; orientation=orientation, frac=overlap_frac
            )

            local_shift = find_pairwise_offset(subA, subB)

            new_ox = ox_neighbor + offsetA_inA[2] + local_shift[1] - offsetB_inB[2]
            new_oy = oy_neighbor + offsetA_inA[1] + local_shift[2] - offsetB_inB[1]

            dist_shift = sqrt((new_ox - ox_guess)^2 + (new_oy - oy_guess)^2)

            if dist_shift <= max_allowed_shift
                push!(local_shifts_x, new_ox)
                push!(local_shifts_y, new_oy)
            end
        end

        if !isempty(local_shifts_x)
            mean_ox = mean(local_shifts_x)
            mean_oy = mean(local_shifts_y)
            dist_shift_final = sqrt((mean_ox - old_ox)^2 + (mean_oy - old_oy)^2)

            if dist_shift_final <= max_allowed_shift
                println("Refined tile (r=$r,c=$c) offset from ",
                    "($old_ox, $old_oy) to ($mean_ox, $mean_oy)")
                refined[tile] = (mean_ox, mean_oy)
            else
                println("**Reverting** tile (r=$r,c=$c): final shift is too big.")
            end
        end
    end

    return refined
end

# ---------------------- align_tiles_reference_timepoint ----------------------
