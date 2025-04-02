# --- Align Functions ---
function tile_using_positions(image_data, bin_factor)
    offsets = Dict{Any,Tuple{Float64,Float64}}()
    for tile in image_data
        px = get(tile, "position_x", 0.0)
        py = get(tile, "position_y", 0.0)
        if px === missing
            px = 0.0
        end
        if py === missing
            py = 0.0
        end
        offsets[tile] = (Float64(px) / bin_factor, Float64(py) / bin_factor)
    end
    return offsets
end


""" Calculates offsets by cross-correlation on one reference plane """
function align_tiles_reference_timepoint(image_data, base_path;
    c_index_ref=1, z_index_ref=1, time_index_ref=1, # Specify which C/Z/T plane to use for alignment
    bin_factor=1, overlap_frac=0.1, max_allowed_shift=100.0)

    println("Aligning based on C=$c_index_ref, Z=$z_index_ref, T=$time_index_ref...")

    # Identify tiles that have the required channel file link
    valid_tiles = [img for img in image_data if
                   haskey(img, "channel_files") &&
                   haskey(img["channel_files"], c_index_ref - 1) && # Check 0-based channel key
                   !ismissing(img["channel_files"][c_index_ref-1]) &&
                   haskey(img, "grid_row") && !ismissing(img["grid_row"]) &&
                   haskey(img, "grid_col") && !ismissing(img["grid_col"])]

    if isempty(valid_tiles)
        println("Error: No valid tiles found for reference alignment C=$c_index_ref, Z=$z_index_ref, T=$time_index_ref.")
        return Dict{Any,Tuple{Float64,Float64}}() # Return empty Z-independent structure
    end

    sort!(valid_tiles, by=x -> (x["grid_row"], x["grid_col"]))

    initial_offsets = Dict{Any,Tuple{Float64,Float64}}() # Z-independent structure
    tile_images_ref = Dict{Any,Matrix{Float32}}() # Store the Float32 binned reference images

    println("Loading reference plane images...")
    # Load & bin the specific reference C/Z/T plane for each valid tile
    for tile in valid_tiles
        c_key = c_index_ref - 1
        fname = tile["channel_files"][c_key]
        path = get_file_path(base_path, fname)
        try
            # load_plane_and_bin handles file finding, loading, Z indexing, binning
            tile_images_ref[tile] = load_plane_and_bin(path, c_index_ref, z_index_ref, time_index_ref, b_factor=bin_factor)
        catch e
            println("Error loading reference plane for tile $(tile["id"]): $e. Skipping tile.")
            # Remove tile from valid_tiles if loading fails? Or handle missing entry later.
            # For now, tile_images_ref will be missing this key.
        end
    end

    # Remove tiles where loading failed from further processing
    filter!(t -> haskey(tile_images_ref, t), valid_tiles)
    if isempty(valid_tiles)
        println("Error: Loading failed for all reference planes.")
        return Dict{Any,Tuple{Float64,Float64}}()
    end

    println("Calculating initial naive offsets...")
    # Naive row/col offset guess based on grid and reference image size
    for tile in valid_tiles
        im = tile_images_ref[tile]
        (h, w) = size(im)
        r = tile["grid_row"]
        c = tile["grid_col"]
        # Use non-integer offsets initially
        ox = Float64(c - 1) * (w * (1.0 - overlap_frac)) # Assume first tile (c=1) is at ox=0
        oy = Float64(r - 1) * (h * (1.0 - overlap_frac)) # Assume first tile (r=1) is at oy=0
        initial_offsets[tile] = (ox, oy)
    end

    println("Refining offsets via local cross-correlation...")
    # Local multi-neighbor refinement using loaded reference images
    # refine_offsets expects Dict{Any, Tuple{Float64,Float64}} and returns the same structure
    refined_offsets = refine_offsets(valid_tiles, tile_images_ref, initial_offsets;
        overlap_frac=overlap_frac,
        max_allowed_shift=max_allowed_shift)

    # refined_offsets now contains the final Z-independent offsets for each tile
    println("Alignment refinement complete.")
    return refined_offsets # Return Dict{ Tile -> (OffsetX, OffsetY) }
end

function overlap_subimages(imA::AbstractMatrix{T}, imB::AbstractMatrix{T}; orientation::Symbol, frac::Float64=0.5) where {T<:Real}
    hA, wA = size(imA)
    hB, wB = size(imB)
    if orientation == :horizontal
        overlap_wA = clamp(round(Int, frac * wA), 1, wA)
        overlap_wB = clamp(round(Int, frac * wB), 1, wB)
        subA = view(imA, :, wA-overlap_wA+1:wA)
        subB = view(imB, :, 1:overlap_wB)
        offsetA = (0, wA - overlap_wA)
        offsetB = (0, 0)
        return subA, subB, offsetA, offsetB
    elseif orientation == :vertical
        overlap_hA = clamp(round(Int, frac * hA), 1, hA)
        overlap_hB = clamp(round(Int, frac * hB), 1, hB)
        subA = view(imA, hA-overlap_hA+1:hA, :)
        subB = view(imB, 1:overlap_hB, :)
        offsetA = (hA - overlap_hA, 0)
        offsetB = (0, 0)
        return subA, subB, offsetA, offsetB
    else
        error("Unknown orientation $orientation")
    end
end

function find_pairwise_offset(imA::AbstractMatrix{T}, imB::AbstractMatrix{T}) where {T<:Real}
    A = Float32.(imA)
    B = Float32.(imB) # Ensure Float32 for FFT
    # Pad smaller image to match larger image size for FFT
    hA, wA = size(A)
    hB, wB = size(B)
    maxH = max(hA, hB)
    maxW = max(wA, wB)
    Ap = A
    Bp = B # Start assuming no padding needed
    if size(A) != (maxH, maxW)
        Ap = PaddedView(0, A, (maxH, maxW))
    end
    if size(B) != (maxH, maxW)
        Bp = PaddedView(0, B, (maxH, maxW))
    end

    FA = fft(Ap)
    FB = fft(Bp)
    R = FA .* conj(FB)
    R ./= abs.(R) .+ 1.0f-6 # Add epsilon for stability
    c = real.(ifft(R))

    peak_idx = argmax(c)
    peak_y, peak_x = Tuple(peak_idx) # Cartesian index unpacks directly
    Ny, Nx = size(c)
    shift_y = peak_y - 1
    shift_x = peak_x - 1 # 0-based shift

    # FFT shift wrapping
    if shift_x > Nx ÷ 2
        shift_x -= Nx
    end
    if shift_y > Ny ÷ 2
        shift_y -= Ny
    end

    return (shift_x, shift_y) # Return (x_shift, y_shift)
end

function refine_offsets(valid_tiles, tile_images::Dict{Any,Matrix{Float32}}, initial_offsets::Dict{Any,Tuple{Float64,Float64}}; overlap_frac::Float64, max_allowed_shift::Float64)
    # Returns Dict{ Any => Tuple{Float64, Float64} }

    refined = copy(initial_offsets)
    tile_map = Dict((t["grid_row"], t["grid_col"]) => t for t in valid_tiles) # Map grid pos to tile object

    for tile in valid_tiles
        if !haskey(refined, tile)
            continue
        end # Skip if no initial offset
        r = tile["grid_row"]
        c = tile["grid_col"]
        (old_ox, old_oy) = refined[tile]
        neighbor_coords = [(r, c - 1), (r, c + 1), (r - 1, c), (r + 1, c)] # Grid coords
        local_shifts_x = Float64[]
        local_shifts_y = Float64[]

        for (nr, nc) in neighbor_coords
            neighbor_tile = get(tile_map, (nr, nc), nothing) # Find neighbor tile object
            if neighbor_tile === nothing || !haskey(refined, neighbor_tile) || !haskey(tile_images, tile) || !haskey(tile_images, neighbor_tile)
                continue # Skip if neighbor doesn't exist, has no offset, or image missing
            end

            # Determine orientation based on grid difference
            orientation = (nr == r) ? :horizontal : :vertical
            # Determine which image is A (left/top) and B (right/bottom)
            imA, imB, tileA, tileB = (c > nc || r > nr) ? (tile_images[neighbor_tile], tile_images[tile], neighbor_tile, tile) : (tile_images[tile], tile_images[neighbor_tile], tile, neighbor_tile)

            (oxA_guess, oyA_guess) = refined[tileA]
            (oxB_guess, oyB_guess) = refined[tileB]

            subA, subB, offsetA_inA, offsetB_inB = overlap_subimages(imA, imB; orientation=orientation, frac=overlap_frac)
            # Check if subimages are valid before proceeding
            if isempty(subA) || isempty(subB)
                println("Warn: Empty subimage for ($r,$c) vs ($nr,$nc)")
                continue
            end

            shift_x_local, shift_y_local = find_pairwise_offset(subA, subB) # Returns (x_shift, y_shift)

            # Calculate the implied absolute position of tile B based on tile A's position and the calculated shift
            # new_pos = posA + offsetA_inA + local_shift - offsetB_inB
            new_ox_B = oxA_guess + offsetA_inA[2] + shift_x_local - offsetB_inB[2]
            new_oy_B = oyA_guess + offsetA_inA[1] + shift_y_local - offsetB_inB[1]

            # Add the calculated absolute position of the *current tile* based on this neighbor pair
            current_tile_is_B = (tile == tileB)
            if current_tile_is_B
                dist_shift = sqrt((new_ox_B - oxB_guess)^2 + (new_oy_B - oyB_guess)^2)
                if dist_shift <= max_allowed_shift
                    push!(local_shifts_x, new_ox_B)
                    push!(local_shifts_y, new_oy_B)
                else
                    println("Warn: Shift too large ($dist_shift > $max_allowed_shift) for ($r,$c) vs ($nr,$nc)")
                end
            else # Current tile is A, we need to infer its position from B's new position
                # posA = posB - offsetA_inA - local_shift + offsetB_inB
                new_ox_A = new_ox_B - offsetA_inA[2] - shift_x_local + offsetB_inB[2]
                new_oy_A = new_oy_B - offsetA_inA[1] - shift_y_local + offsetB_inB[1]
                dist_shift = sqrt((new_ox_A - oxA_guess)^2 + (new_oy_A - oyA_guess)^2)
                if dist_shift <= max_allowed_shift
                    push!(local_shifts_x, new_ox_A)
                    push!(local_shifts_y, new_oy_A)
                else
                    println("Warn: Shift too large ($dist_shift > $max_allowed_shift) for ($r,$c) vs ($nr,$nc)")
                end
            end
        end # End loop over neighbors

        if !isempty(local_shifts_x)
            mean_ox = mean(local_shifts_x)
            mean_oy = mean(local_shifts_y)
            dist_shift_final = sqrt((mean_ox - old_ox)^2 + (mean_oy - old_oy)^2)

            if dist_shift_final <= max_allowed_shift
                # println("Refined ($r,$c) offset from ($old_ox, $old_oy) -> ($mean_ox, $mean_oy)")
                refined[tile] = (mean_ox, mean_oy)
                # else: Keep original offset if refinement suggests too large a shift
                #    println("Warn: ($r,$c) final shift $dist_shift_final too large. Keeping old offset.")
            end
        end
    end # End loop over tiles

    return refined
end