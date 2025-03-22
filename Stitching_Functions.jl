using EzXML
using FileIO

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
                row_val = parse(Int, m.captures[1])
                col_val = parse(Int, m.captures[2])
                img_info["grid_row"] = row_val
                img_info["grid_col"] = col_val
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
        end

        push!(image_data, img_info)
    end

    return image_data
end

function get_actual_file_path(base_dir::String, file_name::String)
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


# ---------------------- bin2 ----------------------
function bin2(img::Matrix{UInt16})
    h, w = size(img)
    h2 = div(h, 2)
    w2 = div(w, 2)
    out = Matrix{UInt16}(undef, h2, w2)

    @inbounds for j in 1:h2
        for i in 1:w2
            v = UInt32(img[2j-1, 2i-1]) +
                UInt32(img[2j-1, 2i])   +
                UInt32(img[2j,   2i-1]) +
                UInt32(img[2j,   2i])
            v >>= 2
            v = min(v, 0xFFFF)
            out[j, i] = UInt16(v)
        end
    end
    return out
end



# ---------------------- load_and_bin ----------------------
function load_and_bin(path::String; time_index::Int=1, channel_index::Int=1)
    full_img = TiffImages.load(path)
    nd = ndims(full_img)

    plane = nothing
    if nd == 2
        plane = full_img
    elseif nd == 3
        @assert time_index ≤ size(full_img, 3) "time_index out of range"
        plane = full_img[:, :, time_index]
    elseif nd == 4
        @assert channel_index ≤ size(full_img, 3) "channel_index out of range"
        @assert time_index   ≤ size(full_img, 4) "time_index out of range"
        plane = full_img[:, :, channel_index, time_index]
    else
        throw(ArgumentError("Unsupported dimension $(ndims(full_img)) for TIF at $path"))
    end

    h, w = size(plane)
    tmp = Matrix{UInt16}(undef, h, w)

    @inbounds for row in 1:h
        for col in 1:w
            val = plane[row, col]
            floatval = Float32(val) * 65535f0
            clamped  = clamp(round(UInt32, floatval), 0, 0xFFFF)
            tmp[row, col] = UInt16(clamped)
        end
    end

    return bin2(tmp)
end