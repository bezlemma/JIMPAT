function parse_ome(file_path::String)
    doc = readxml(file_path)
    root_el = root(doc)

    df = DataFrame(
        pos_id     = Int[],
        c          = Int[],
        z          = Int[],
        t          = Int[],
        filename   = String[],
        x_position = Float64[],
        y_position = Float64[]
    )

    ns = Dict("ome" => "http://www.openmicroscopy.org/Schemas/OME/2016-06")
    image_nodes = findall("//ome:Image", root_el, ns)
    
    # Pre-scan grid dimensions (if grid notation is used) by finding the max column value.
    grid_width = nothing
    for image_elem in image_nodes
        stage_label = findfirst(".//ome:StageLabel", image_elem, ns)
        if stage_label !== nothing && haskey(stage_label, "Name")
            label_text = stage_label["Name"]
            if occursin("Row", label_text)
                m = match(r"Row(\d+)_Col(\d+)", label_text)
                if m !== nothing
                    col = parse(Int, m.captures[2])
                    grid_width = isnothing(grid_width) ? (col + 1) : max(grid_width, col + 1)
                end
            end
        end
    end

    for image_elem in image_nodes
        # Determine pos_id based on StageLabel.
        pos_id = 0
        stage_label = findfirst(".//ome:StageLabel", image_elem, ns)
        if stage_label !== nothing && haskey(stage_label, "Name")
            label_text = stage_label["Name"]
            if occursin("Number", label_text)
                m = match(r"Number(\d+)_sg", label_text)
                if m !== nothing
                    pos_id = parse(Int, m.captures[1])
                end
            elseif occursin("Row", label_text)
                m = match(r"Row(\d+)_Col(\d+)", label_text)
                if m !== nothing
                    row = parse(Int, m.captures[1])
                    col = parse(Int, m.captures[2])
                    local_gw = isnothing(grid_width) ? (col + 1) : grid_width
                    pos_id = row * local_gw + col + 1
                end
            end
        end

        pixels_elem = findfirst(".//ome:Pixels", image_elem, ns)
        size_z = parse(Int, pixels_elem["SizeZ"])
        size_c = parse(Int, pixels_elem["SizeC"])
        size_t = parse(Int, pixels_elem["SizeT"])
        px_size_x = haskey(pixels_elem, "PhysicalSizeX") ? parse(Float64, pixels_elem["PhysicalSizeX"]) : 1.0
        px_size_y = haskey(pixels_elem, "PhysicalSizeY") ? parse(Float64, pixels_elem["PhysicalSizeY"]) : 1.0

        # Build a dictionary (c, t, z) => filename
        file_map = Dict{Tuple{Int,Int,Int}, String}()

        tiff_data_nodes = findall(".//ome:TiffData", pixels_elem, ns)
        for td in tiff_data_nodes
            uuid_elem = findfirst(".//ome:UUID", td, ns)
            if uuid_elem === nothing || !haskey(uuid_elem, "FileName")
                continue
            end

            fname = uuid_elem["FileName"]

            # If an attribute is absent, default to 0
            fc = haskey(td, "FirstC") ? parse(Int, td["FirstC"]) : 0
            ft = haskey(td, "FirstT") ? parse(Int, td["FirstT"]) : 0
            fz = haskey(td, "FirstZ") ? parse(Int, td["FirstZ"]) : 0

            # Store single mapping: (fc, ft, fz) => fname
            file_map[(fc, ft, fz)] = fname
        end

        # Now handle each <Plane>, look up (c0,t0,z0) in the dictionary
        plane_nodes = findall(".//ome:Plane", pixels_elem, ns)
        for plane_node in plane_nodes
            c0 = parse(Int, haskey(plane_node, "TheC") ? plane_node["TheC"] : "0")
            t0 = parse(Int, haskey(plane_node, "TheT") ? plane_node["TheT"] : "0")
            z0 = parse(Int, haskey(plane_node, "TheZ") ? plane_node["TheZ"] : "0")

            c_ = c0 + 1
            t_ = t0 + 1
            z_ = z0 + 1

            px = haskey(plane_node, "PositionX") ? parse(Float64, plane_node["PositionX"]) : 0.0
            py = haskey(plane_node, "PositionY") ? parse(Float64, plane_node["PositionY"]) : 0.0
            x_ = px_size_x != 0.0 ? px / px_size_x : 0.0
            y_ = px_size_y != 0.0 ? py / px_size_y : 0.0

            the_file = get(file_map, (c0, t0, z0), nothing)
            if the_file === nothing
                error("No filename found for c=$c0, t=$t0, z=$z0 at pos_id=$pos_id.")
            end

            push!(df, (pos_id, c_, z_, t_, the_file, x_, y_))
        end
    end

    return df
end





function parse_ome_metadata(file_path::String)
    doc = readxml(file_path)
    root_el = root(doc)
    ns = Dict("ome" => "http://www.openmicroscopy.org/Schemas/OME/2016-06")

    image_elem = findfirst("//ome:Image", root_el, ns)
    pixels_elem = findfirst(".//ome:Pixels", image_elem, ns)

    meta_info = Dict{String, Any}()

    meta_info["SizeZ"] = parse(Int, pixels_elem["SizeZ"])
    meta_info["SizeC"] = parse(Int, pixels_elem["SizeC"])
    meta_info["SizeT"] = parse(Int, pixels_elem["SizeT"])
    meta_info["SizeX"] = parse(Int, pixels_elem["SizeX"])
    meta_info["SizeY"] = parse(Int, pixels_elem["SizeY"])
    meta_info["PhysicalSizeX"] = haskey(pixels_elem, "PhysicalSizeX") ? parse(Float64, pixels_elem["PhysicalSizeX"]) : 1.0
    meta_info["PhysicalSizeY"] = haskey(pixels_elem, "PhysicalSizeY") ? parse(Float64, pixels_elem["PhysicalSizeY"]) : 1.0
    meta_info["PhysicalSizeZ"] = haskey(pixels_elem, "PhysicalSizeZ") ? parse(Float64, pixels_elem["PhysicalSizeZ"]) : 1.0
    meta_info["DimensionOrder"] = haskey(pixels_elem, "DimensionOrder") ? pixels_elem["DimensionOrder"] : "TCZYX"
    meta_info["Type"] = haskey(pixels_elem, "Type") ? pixels_elem["Type"] : "uint16"
    meta_info["BigEndian"] = haskey(pixels_elem, "BigEndian") ? pixels_elem["BigEndian"] : "false"

    return meta_info
end



function load_image(image_data::DataFrame, base_path::String, position_id, c_index::Int, z_index::Int, t_index::Int)
    row = filter(r -> r.pos_id == position_id && r.c == c_index && r.z == z_index && r.t == t_index, image_data)
    if nrow(row) == 0
        error("No matching entry found for pos_id=$(position_id), c=$(c_index), z=$(z_index), t=$(t_index)")
    end
    filename = row[1, :filename]
    path = isabspath(filename) ? filename : joinpath(base_path, filename)
    raw_tiff_stack = TiffImages.load(path; verbose=false)
    if ndims(raw_tiff_stack) == 3
        raw_plane = raw_tiff_stack[:, :, z_index]
    elseif ndims(raw_tiff_stack) == 2
        raw_plane = raw_tiff_stack
    else
        error("Loaded stack from '$path' has unexpected dimensions: $(ndims(raw_tiff_stack))")
    end
    return raw_plane
end


function bin2(img)
    h, w = size(img)
    h2 = div(h, 2)
    w2 = div(w, 2)
    out = similar(img, h2, w2)
    for j in 1:h2, i in 1:w2
        s = img[2j-1, 2i-1] + img[2j-1, 2i] + img[2j, 2i-1] + img[2j, 2i]
        out[j, i] = s / 4
    end
    return out
end

function bin4(img)
    h, w = size(img)
    h4 = div(h, 4)
    w4 = div(w, 4)
    out = similar(img, h4, w4)
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