using TiffImages        # for loading & saving TIF volumes
using ImageIO           # alternative load/save
using ImageMorphology   # for 3D dilation, etc.
using ImageSegmentation # for connected components
using Statistics        # for mean, cov, etc. (used in PCA step)
using LinearAlgebra     # for eigen decomposition
using Colors

# ---------------------------
## 1) Define paths and parameters
# ---------------------------
PATH         = raw"D:\MyFirstPaperData\10x_Fixed_NoCulture_MAICO\batch1_done"
OGDATA_PATH  = "try2_10xhighna_2branch_hourly_2p56z.tif"
BINARY_PATH  = "MASK_ECad.tif"
CHANNELS     = 4

# ---------------------------
## 2) Load the data
# ---------------------------
BINARY_DATA = TiffImages.load(joinpath(PATH, BINARY_PATH))  # Might be UInt8/UInt16
OG_DATA     = TiffImages.load(joinpath(PATH, OGDATA_PATH))   # Multi-channel TIF

# Convert BINARY_DATA to Bool if needed (Ilastik often produces 0..255)
BINARY_DATA = BINARY_DATA .> 0

# Separate channels in OG_DATA if channels are interleaved in the 3rd dimension
# (Assumes shape is (height, width, colorSlices). Adjust if needed.)
OG_DATA_CH1 = OG_DATA[:, :, 1:CHANNELS:end]
OG_DATA_CH2 = CHANNELS > 1 ? OG_DATA[:, :, 2:CHANNELS:end] : nothing
OG_DATA_CH3 = CHANNELS > 2 ? OG_DATA[:, :, 3:CHANNELS:end] : nothing
OG_DATA_CH4 = CHANNELS > 3 ? OG_DATA[:, :, 4:CHANNELS:end] : nothing

# ---------------------------
## 3) Morphological Dilation & Fill Holes (3D)
# ---------------------------
# 3D structuring element: 3×3×3 "sphere" approximation
kernel = ones(Bool, 3, 3, 3)

# Dilation
BINARY_DATA = ImageMorphology.dilate(BINARY_DATA, kernel)

# Fill holes in 3D by labeling the complement & inverting
function fill_3d_holes(vol)
    invVol = .!vol
    # Fixed: use structuring element instead of connectivity parameter
    se = ImageMorphology.strel(ones(Bool, 3, 3, 3))
    labeled = ImageMorphology.label_components(invVol, se)
    
    outside_label = labeled[1,1,1]
    filled = copy(vol)
    filled[labeled .!= outside_label] .= true
    return filled
end
BINARY_DATA = fill_3d_holes(BINARY_DATA)

# ---------------------------
## 4) Keep Only Largest Connected Component
# ---------------------------
# Fixed: use structuring element for connectivity
se = ImageMorphology.strel(ones(Bool, 3, 3, 3))
labeled = ImageMorphology.label_components(BINARY_DATA, se)

# Count voxels for each label
label_counts = Dict{Int, Int}()
for i in labeled
    if i != 0  # Skip background (0)
        label_counts[i] = get(label_counts, i, 0) + 1
    end
end

# Find label with most voxels
function find_biggest_label(label_counts)
    biggest_label = 0
    max_count = 0
    for (label, count) in label_counts
        if count > max_count
            max_count = count
            biggest_label = label
        end
    end
    return biggest_label
end

biggest_label = find_biggest_label(label_counts)
# Keep only the largest component
BINARY_DATA .= (labeled .== biggest_label)

# ---------------------------
## 5) Multiply the Final Binary Mask by the Original Data
# ---------------------------
BINARY_DATA_U16 = convert.(UInt16, BINARY_DATA)
NEW_CH1 = BINARY_DATA_U16 .* OG_DATA_CH1
NEW_CH2 = CHANNELS > 1 ? BINARY_DATA_U16 .* OG_DATA_CH2 : nothing
NEW_CH3 = CHANNELS > 2 ? BINARY_DATA_U16 .* OG_DATA_CH3 : nothing
NEW_CH4 = CHANNELS > 3 ? BINARY_DATA_U16 .* OG_DATA_CH4 : nothing

# ---------------------------
## 6) Approximate Volume, Surface Area, and Principal Axis
# ---------------------------
MAG       = 10.0
BIN       = 1.0
Z_STEP    = 6.25
um2Pixel  = 6.25   # microscope scale
pixel2len = um2Pixel * BIN / MAG
pixel2area   = pixel2len^2
voxel2volume = pixel2len^2 * Z_STEP

# Volume
voxelcount = count(BINARY_DATA)

# Surface area (approximation by counting boundary faces)
function voxel_surface_count(vol)
    (nx, ny, nz) = size(vol)
    s = 0
    # check faces along x
    for x in 1:(nx-1), y in 1:ny, z in 1:nz
        s += Int(vol[x,y,z] != vol[x+1,y,z])
    end
    # check faces along y
    for x in 1:nx, y in 1:(ny-1), z in 1:nz
        s += Int(vol[x,y,z] != vol[x,y+1,z])
    end
    # check faces along z
    for x in 1:nx, y in 1:ny, z in 1:(nz-1)
        s += Int(vol[x,y,z] != vol[x,y,z+1])
    end
    return s
end
surfcount = voxel_surface_count(BINARY_DATA)

# Principal axis length (3D PCA on voxel coords)
function principal_axis_length(vol)
    # Fixed: collect coordinates more efficiently
    coords = Array{Float64}(undef, count(vol), 3)
    idx = 1
    for I in CartesianIndices(vol)
        if vol[I]
            coords[idx, 1] = I[1]
            coords[idx, 2] = I[2]
            coords[idx, 3] = I[3]
            idx += 1
        end
    end
    
    cmean = mean(coords, dims=1)
    ccentered = coords .- cmean
    Σ = cov(ccentered)
    evals, evecs = eigen(Σ)
    
    # Find index of maximum eigenvalue
    max_idx = argmax(evals)
    mainvec = evecs[:, max_idx]
    
    projs = ccentered * mainvec
    return maximum(projs) - minimum(projs)
end

paxis_len_voxels = principal_axis_length(BINARY_DATA)

vol_mm3      = voxel2volume * voxelcount / (1000^3)
surf_mm2     = pixel2area   * surfcount   / (1000^2)
paxis_len_mm = pixel2len    * paxis_len_voxels / 1000

println("Volume is $vol_mm3 mm^3")
println("Surface Area is $surf_mm2 mm^2")
println("Length is $paxis_len_mm mm")

# ---------------------------
## 7) Save the Final Mask and Channel Data
# ---------------------------
function save_tiff_stack(filename::String, vol::BitArray{3})
    # Convert BitArray to UInt8 (0 and 1)
    vol_uint8 = UInt8.(vol)
    # Convert each pixel to Gray{N0f8} (0 becomes black, 1 becomes white)
    vol_gray = map(x -> Gray{N0f8}(x == 1 ? 1.0 : 0.0), vol_uint8)
    # Save the 3D array as a multi-page TIFF
    TiffImages.save(filename, vol_gray)
end
PATH = "./"
save_tiff_stack(joinpath(PATH, "newseg.tif"), BINARY_DATA)

if NEW_CH1 !== nothing
    save_tiff_stack(joinpath(PATH, "newdata_CH1.tif"), NEW_CH1)
end
if NEW_CH2 !== nothing
    save_tiff_stack(joinpath(PATH, "newdata_CH2.tif"), NEW_CH2)
end
if NEW_CH3 !== nothing
    save_tiff_stack(joinpath(PATH, "newdata_CH3.tif"), NEW_CH3)
end
if NEW_CH4 !== nothing
    save_tiff_stack(joinpath(PATH, "newdata_CH4.tif"), NEW_CH4)
end

println("Processing complete. Files saved to $PATH")

## View, Plot 3D

using GLMakie

figure = Figure()
axis3 = Axis3(figure[1, 1],  aspect = (1, 1, 1),)

colormap = to_colormap(:plasma)
colormap[1] = RGBAf(0,0,0,0)

volume!(BINARY_DATA;
    algorithm   = :absorption,
    colormap=colormap
)

display(figure)