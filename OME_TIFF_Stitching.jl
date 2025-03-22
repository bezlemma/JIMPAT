using FileIO, FFTW, Base.Threads, LinearAlgebra, TiffImages
using Statistics
using SparseArrays  # Needed for the global solver
using OrderedCollections

include("./Stitching_Functions.jl")

function overlap_subimages(imA::Matrix{UInt16}, imB::Matrix{UInt16};
                           orientation::Symbol, frac::Float64=0.5)
    hA, wA = size(imA)
    hB, wB = size(imB)

    if orientation == :horizontal
        overlap_widthA = round(Int, frac * wA)
        overlap_widthB = round(Int, frac * wB)
        subA = imA[:, wA - overlap_widthA + 1 : wA]
        subB = imB[:, 1 : overlap_widthB]
        offsetA = (0, wA - overlap_widthA)
        offsetB = (0, 0)
        return subA, subB, offsetA, offsetB

    elseif orientation == :vertical
        overlap_heightA = round(Int, frac * hA)
        overlap_heightB = round(Int, frac * hB)
        subA = imA[hA - overlap_heightA + 1 : hA, :]
        subB = imB[1 : overlap_heightB, :]
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
    R ./= abs.(R) .+ 1f0
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
        neighbor_positions = [(r, c-1), (r, c+1), (r-1, c), (r+1, c)]
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
            (ox_guess,   oy_guess)     = refined[tile]

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
function align_tiles_reference_timepoint(image_data, base_path;
                                         time_index=0, overlap_frac=0.1,
                                         max_allowed_shift=100.0)
    valid_tiles = [img for img in image_data if
        img["grid_row"] !== missing &&
        img["grid_col"] !== missing &&
        haskey(img["tiff_files"], time_index)
    ]
    sort!(valid_tiles, by = x->(x["grid_row"], x["grid_col"]))

    offsets = Dict{Any, Tuple{Float64, Float64}}()
    tile_images = Dict{Any, Matrix{UInt16}}()

    # load & bin each tile
    for tile in valid_tiles
        fname = tile["tiff_files"][time_index]
        path = get_actual_file_path(base_path, fname)
        tile_images[tile] = load_and_bin(path; time_index = time_index + 1)
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
        sort!(row_tiles, by = x->x["grid_col"])
        row_imgs = [ tile_images[t] for t in row_tiles ]
        push!(row_blocks, hcat(row_imgs...))
    end
    rowcol_image = vcat(row_blocks...)

    return refined_offsets, rowcol_image
end

# ---------------------- Feathered blending version ----------------------
function build_feather_mask(h::Int, w::Int; feather_size::Float64=30.0)::Matrix{Float32}
    mask = Matrix{Float32}(undef, h, w)
    for j in 1:h
        dist_top    = j - 1
        dist_bottom = h - j
        for i in 1:w
            dist_left  = i - 1
            dist_right = w - i
            min_edge_dist = min(dist_top, dist_bottom, dist_left, dist_right)
            if min_edge_dist >= feather_size
                mask[j, i] = 1.0f0
            else
                # transitions linearly from 1.0 to 0.0 over 'feather_size' near the edges
                mask[j, i] = Float32(min_edge_dist / feather_size)
            end
        end
    end
    return mask
end

function stitch_timepoint(image_data, offsets, time_index::Int, base_path::String;
                          feather_size::Float64=30.0)
    relevant_tiles = [img for img in image_data if haskey(offsets, img)]
    if isempty(relevant_tiles)
        return fill(UInt16(0), 1, 1)
    end

    tile_images = Dict{Any, Matrix{UInt16}}()
    xvals = Float64[]
    yvals = Float64[]

    # load & bin each tile
    for tile in relevant_tiles
        if haskey(tile["tiff_files"], time_index)
            fname = tile["tiff_files"][time_index]
            path = get_actual_file_path(base_path, fname)
            binned_img = load_and_bin(path; time_index = time_index + 1)
            tile_images[tile] = binned_img

            (ox, oy) = offsets[tile]
            w = size(binned_img, 2)
            h = size(binned_img, 1)
            push!(xvals, ox, ox + w - 1)
            push!(yvals, oy, oy + h - 1)
        end
    end

    if isempty(xvals)
        return fill(UInt16(0), 1, 1)
    end

    min_x = floor(Int, minimum(xvals))
    max_x = ceil(Int, maximum(xvals))
    min_y = floor(Int, minimum(yvals))
    max_y = ceil(Int, maximum(yvals))
    width  = max_x - min_x + 1
    height = max_y - min_y + 1

    # We'll accumulate weighted sums here and then divide by total weight
    accum_image  = fill(0.0f0, height, width)
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

        h_tile = size(im, 1)
        w_tile = size(im, 2)

        # build a feather mask for this tile
        mask = build_feather_mask(h_tile, w_tile; feather_size=feather_size)

        @inbounds for j in 1:h_tile
            for i in 1:w_tile
                wgt = mask[j, i]
                if wgt > 0
                    accum_image[oy_s + j, ox_s + i]  += float(im[j, i]) * wgt
                    accum_weight[oy_s + j, ox_s + i] += wgt
                end
            end
        end
    end

    # Convert accumulated floats back to UInt16
    canvas = fill(UInt16(0), height, width)
    @inbounds for y in 1:height
        for x in 1:width
            wgt = accum_weight[y, x]
            if wgt > 0
                val = accum_image[y, x] / wgt
                val_clamped = clamp(round(Int, val), 0, 65535)
                canvas[y, x] = UInt16(val_clamped)
            end
        end
    end

    return canvas
end

# ---------------------- Main script usage example ----------------------
companion_file = raw"C:\Users\uComp\Documents\LinneaData\20250226\Sample5_1\20250226_ll_sample5_60x_10msexp_1-2bint_gfp-bypass_60sint1.companion.ome"
base_path      = raw"C:\Users\uComp\Documents\LinneaData\20250226\Sample5_1"
output_file    = raw"C:\Users\uComp\Documents\LinneaData\20250226\Sample5_1\output_Julia_binned.tif"

image_data = parse_ome_companion(companion_file)
max_t = maximum(img["size_t"] for img in image_data)  # largest 0-based index
ref_time = max_t - 1
println("Aligning tiles based on final timepoint (index $ref_time).")

(offsets, rowcol_image) = align_tiles_reference_timepoint(image_data, base_path;
    time_index       = ref_time,
    overlap_frac     = 0.1,
    max_allowed_shift= 200.0
)
println("Done computing & refining offsets at final timepoint.")

println("Stitching reference time point (index $ref_time) first ...")
mosaic_ref = stitch_timepoint(image_data, offsets, ref_time, base_path; feather_size=30.0)
(h0, w0) = size(mosaic_ref)

stitched_series = Array{UInt16,3}(undef, h0, w0, max_t)
stitched_series[:,:,ref_time+1] = mosaic_ref

@threads for t in 40:50  # you can adjust this range as needed
    if t == ref_time
        continue
    end
    println("Stitching time point $(t+1) of $max_t ...")
    mosaic_t = stitch_timepoint(image_data, offsets, t, base_path; feather_size=30.0)
    if size(mosaic_t) != (h0, w0)
        error("Time $t mosaic dimension != reference mosaic. Possibly missing tiles or mismatch.")
    end
    stitched_series[:,:,t+1] = mosaic_t
end

save(output_file, stitched_series)
println("Saved stitched binned series (UInt16) to: $output_file")
println("Done!")
