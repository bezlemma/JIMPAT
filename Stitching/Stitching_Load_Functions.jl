using EzXML, FileIO

function parse_ome_companion(file_path::String)
    doc = readxml(file_path)
    root_node = root(doc)
    ome_ns = namespace(root_node)
    ns_dict = isempty(ome_ns) ? Dict{String,String}() : ["ns" => ome_ns]
    image_nodes = findall("//ns:Image", root_node, ns_dict)
    image_data = Dict{String,Any}[]

    for image_elem in image_nodes
        img_info = Dict{String,Any}()
        img_info["id"] = haskey(image_elem, "ID") ? image_elem["ID"] : "UnknownID_$(length(image_data))"
        img_info["name"] = haskey(image_elem, "Name") ? image_elem["Name"] : "UnknownName_$(length(image_data))"

        stage_label_elem = findfirst("ns:StageLabel", image_elem, ns_dict)
        label_text = ""
        if stage_label_elem !== nothing && haskey(stage_label_elem, "Name")
            label_text = stage_label_elem["Name"]
        end

        m = match(r"(?:Row|R)?(\d+)[_\s,-]*(?:Col|C)?(\d+)", label_text) # Try Row/Col first
        if m !== nothing
            img_info["grid_row"] = parse(Int, m.captures[1])
            img_info["grid_col"] = parse(Int, m.captures[2])
        else
            m_alt = match(r"(\d+):Number(\d+)_sg:(\d+)", label_text) # Try specific pattern from metadata
            if m_alt !== nothing
                img_info["grid_row"] = parse(Int, m_alt.captures[2]) # Use tile index as row
                img_info["grid_col"] = 1 # Assume single column
            else
                println("Warning: Could not parse Row/Col or Index from StageLabel: '$label_text'. Using tile order.")
                img_info["grid_row"] = length(image_data) + 1 # Fallback to tile order
                img_info["grid_col"] = 1
            end
        end

        # --- Pixels Element Parsing ---
        pixels_elem = findfirst("ns:Pixels", image_elem, ns_dict)

        img_info["size_x"] = parse(Int, pixels_elem["SizeX"])
        img_info["size_y"] = parse(Int, pixels_elem["SizeY"])
        img_info["metadata_size_z"] = parse(Int, pixels_elem["SizeZ"])
        img_info["metadata_size_c"] = parse(Int, pixels_elem["SizeC"])
        img_info["metadata_size_t"] = parse(Int, pixels_elem["SizeT"])
        img_info["pixel_type"] = pixels_elem["Type"]

        # Get DimensionOrder, default if missing
        if haskey(pixels_elem, "DimensionOrder")
            img_info["metadata_dimension_order"] = uppercase(pixels_elem["DimensionOrder"])
        else
            println("Warning: Image $(img_info["id"]) - Pixels element missing 'DimensionOrder' attribute. Assuming 'XYZCT'.")
            img_info["metadata_dimension_order"] = "XYZCT"
        end

        px_size_x = haskey(pixels_elem, "PhysicalSizeX") ? parse(Float64, pixels_elem["PhysicalSizeX"]) : 1.0
        px_size_y = haskey(pixels_elem, "PhysicalSizeY") ? parse(Float64, pixels_elem["PhysicalSizeY"]) : 1.0
        img_info["pixel_size_x"] = px_size_x
        img_info["pixel_size_y"] = px_size_y

        # --- TiffData Parsing: Map Channel Index to Filename ---
        img_info["channel_files"] = Dict{Int,String}() # Map: 0-based C index -> filename
        tiff_data_nodes = findall("ns:TiffData", pixels_elem, ns_dict)

        processed_channels = Set{Int}()
        first_tiffdata = true
        for td in tiff_data_nodes
            uuid_elem = findfirst("ns:UUID", td, ns_dict)
            if uuid_elem !== nothing && haskey(uuid_elem, "FileName")
                filename = uuid_elem["FileName"]
                c_str = haskey(td, "FirstC") ? td["FirstC"] : (first_tiffdata ? "0" : missing)
                first_tiffdata = false
                if ismissing(c_str)
                    continue
                end
                c_index_0based = parse(Int, c_str)
                if !haskey(img_info["channel_files"], c_index_0based)
                    img_info["channel_files"][c_index_0based] = filename
                    push!(processed_channels, c_index_0based)
                end

            end
        end
        expected_channels = 0:(img_info["metadata_size_c"]-1)
        if length(processed_channels) != img_info["metadata_size_c"] && !isempty(processed_channels)
            missing_channels = setdiff(expected_channels, keys(img_info["channel_files"]))
            # println("Warning: Image $(img_info["id"]): File links missing for channels: $(collect(missing_channels)).")
        end
        if isempty(img_info["channel_files"])
            println("Error: Image $(img_info["id"]): No valid channel file links found in TiffData.")
        end

        # --- Plane Position Parsing ---
        plane_elem = findfirst("ns:Plane", pixels_elem, ns_dict)
        if plane_elem !== nothing
            posx = haskey(plane_elem, "PositionX") ? parse(Float64, plane_elem["PositionX"]) : 0.0
            posy = haskey(plane_elem, "PositionY") ? parse(Float64, plane_elem["PositionY"]) : 0.0
            img_info["position_x"] = px_size_x != 0 ? (posx / px_size_x) : 0.0
            img_info["position_y"] = px_size_y != 0 ? (posy / px_size_y) : 0.0
        else
            println("Warning: No Plane element for position in Image $(img_info["id"]). Assuming (0,0).")
            img_info["position_x"] = 0.0
            img_info["position_y"] = 0.0
        end


        push!(image_data, img_info)

    end
    return image_data

end

"""
Loads a specific plane from a TIFF file
"""
#TODO: Make this equally happy to load t, z or c 
function load_plane_and_bin(path::String, c_index::Int, z_index::Int, t_index::Int; b_factor::Int=1)

    plane_idx = max(z_index, t_index) #Little hack right now, likely one of these is 0

    raw_tiff_stack = TiffImages.load(path; verbose=false)
    stack_dims = ndims(raw_tiff_stack)
    local raw_plane


    if stack_dims == 3 #we've loaded a stack, and we need to pick a place
        raw_plane = raw_tiff_stack[:, :, plane_idx]
    elseif stack_dims == 2 #we've loaded a file, and that is all there is to it
        raw_plane = raw_tiff_stack
    else
        error("Loaded stack unexpected dims: $(size(raw_tiff_stack)). File: $path")
    end

    plane_float = Float32.(raw_plane)

    if b_factor == 1
        return plane_float
    elseif b_factor == 2
        return bin2(plane_float)
    elseif b_factor == 4
        return bin4(plane_float)
    else
        println("Warn: Unsupported bin factor $b_factor.")
        return plane_float
    end

end

# --- Binning Functions ---
function bin2(img::AbstractMatrix{T}) where {T<:Real}
    h, w = size(img)
    h2 = h / 2;    w2 = w / 2
    out = Matrix{Float32}(undef, h2, w2)
    for j in 1:h2, i in 1:w2
        s = Float32(img[2j-1, 2i-1]) + Float32(img[2j-1, 2i]) + Float32(img[2j, 2i-1]) + Float32(img[2j, 2i])
        out[j, i] = s / 4.0f0
    end
    return out
end

function bin4(img::AbstractMatrix{T}) where {T<:Real}
    h, w = size(img)
    h4 = h / 4;    w4 = w / 4
    out = Matrix{Float32}(undef, h4, w4)
    for j in 1:h4, i in 1:w4
        s = 0.0f0
        for row_offset in 0:3, col_offset in 0:3
            s += Float32(img[4j-3+row_offset, 4i-3+col_offset])
        end
        out[j, i] = s / 16.0f0
    end
    return out
end

# --- File Path Helper ---
function get_file_path(base_dir::String, file_name::String)
    abs_path = isabspath(file_name) ? file_name : joinpath(base_dir, file_name)
    if isfile(abs_path)
        return abs_path
    end
    try
        for f in readdir(base_dir)
            if lowercase(f) == lowercase(basename(file_name))
                found_path = joinpath(base_dir, f)
                if isfile(found_path)
                    return found_path
                end
            end
        end
    catch e
        println("Warn: Cannot read dir $base_dir: $e")
    end
    return abs_path
end