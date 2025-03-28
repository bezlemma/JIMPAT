using EzXML, FileIO

function parse_ome_companion(file_path::String)
    doc = readxml(file_path)
    root_node = root(doc)

    ome_ns = namespace(root_node)
    ns_dict = ["ns" => ome_ns]

    image_nodes = findall("//ns:Image", root_node, ns_dict)
    image_data = Dict{String,Any}[]

    for image_elem in image_nodes
        img_info = Dict{String,Any}()
        img_info["id"]   = image_elem["ID"]
        img_info["name"] = image_elem["Name"]

        stage_label_elem = findfirst("ns:StageLabel", image_elem, ns_dict)
        if stage_label_elem !== nothing
            label_text = stage_label_elem["Name"]
            m = match(r"(?:\d+:)?Row(\d+)_Col(\d+)(?:_sg:\d+)?", label_text)
            if m !== nothing
                img_info["grid_row"] = parse(Int, m.captures[1])
                img_info["grid_col"] = parse(Int, m.captures[2])
            else
                img_info["grid_row"] = missing
                img_info["grid_col"] = missing
            end
        else
            img_info["grid_row"] = missing
            img_info["grid_col"] = missing
        end

        pixels_elem = findfirst("ns:Pixels", image_elem, ns_dict)
        if pixels_elem !== nothing
            img_info["size_x"] = parse(Int, pixels_elem["SizeX"])
            img_info["size_y"] = parse(Int, pixels_elem["SizeY"])
            img_info["size_z"] = parse(Int, pixels_elem["SizeZ"])
            img_info["size_t"] = parse(Int, pixels_elem["SizeT"])
            img_info["size_c"] = parse(Int, pixels_elem["SizeC"])

            px_size_x = haskey(pixels_elem, "PhysicalSizeX") ? parse(Float64, pixels_elem["PhysicalSizeX"]) : 1.0
            px_size_y = haskey(pixels_elem, "PhysicalSizeY") ? parse(Float64, pixels_elem["PhysicalSizeY"]) : 1.0
            img_info["pixel_size_x"] = px_size_x
            img_info["pixel_size_y"] = px_size_y

            tiff_map = Dict{Int,String}()
            for td in findall("ns:TiffData", pixels_elem, ns_dict)
                t_index = haskey(td, "FirstT") ? parse(Int, td["FirstT"]) : 0
                uuid_elem = findfirst("ns:UUID", td, ns_dict)
                if uuid_elem !== nothing
                    tiff_map[t_index] = uuid_elem["FileName"]
                end
            end
            img_info["tiff_files"] = tiff_map

            plane_elem = findfirst("ns:Plane", pixels_elem, ns_dict)
            if plane_elem !== nothing
                posx = haskey(plane_elem, "PositionX") ? parse(Float64, plane_elem["PositionX"]) : 0.0
                posy = haskey(plane_elem, "PositionY") ? parse(Float64, plane_elem["PositionY"]) : 0.0
                img_info["position_x"] = posx / px_size_x
                img_info["position_y"] = posy / px_size_y
            else
                img_info["position_x"] = missing
                img_info["position_y"] = missing
            end
        end

        push!(image_data, img_info)
    end

    return image_data
end


function load_and_bin(path::String; b_factor::Int=1)
    plane = TiffImages.load(path)
    if b_factor == 1
        return plane
    elseif b_factor == 4
        return bin4(plane)
    else
        return bin2(plane)
    end
end

function bin2(img)
    h, w = size(img)
    h2 = div(h, 2)
    w2 = div(w, 2)
    out = Matrix{Float32}(undef, h2, w2)
    @inbounds for j in 1:h2, i in 1:w2
        s = Float32(img[2j-1, 2i-1]) +
            Float32(img[2j-1, 2i])   +
            Float32(img[2j,   2i-1]) +
            Float32(img[2j,   2i])
        out[j, i] = s / 4.0f0
    end
    return out
end

function bin4(img)
    h, w = size(img)
    h4 = div(h, 4)
    w4 = div(w, 4)
    out = Matrix{Float32}(undef, h4, w4)
    @inbounds for j in 1:h4, i in 1:w4
        s = Float32(img[4j-1, 4i-1]) +
            Float32(img[4j-1, 4i])   +
            Float32(img[4j,   4i-1]) +
            Float32(img[4j,   4i])
        out[j, i] = s / 16.0f0
    end
    return out
end

function get_file_path(base_dir::String, file_name::String)
    p = joinpath(base_dir, file_name)
    if isfile(p)
        return p
    end
    for f in readdir(base_dir)
        if lowercase(f) == lowercase(file_name)
            return joinpath(base_dir, f)
        end
    end
    return p
end