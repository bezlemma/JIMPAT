using TiffImages, ImageIO, ImageMorphology
using ImageSegmentation, Statistics, LinearAlgebra
using GLMakie

#parameters
PATH         = raw"D:\MyFirstPaperData\10x_Fixed_NoCulture_MAICO\batch1_done"
OGDATA_PATH  = "try2_10xhighna_2branch_hourly_2p56z.tif"
BINARY_PATH  = "MASK_ECad.tif"
CHANNELS     = 4

## Load the data
BINARY_DATA = TiffImages.load(joinpath(PATH, BINARY_PATH))  # Might be UInt8/UInt16
OG_DATA     = TiffImages.load(joinpath(PATH, OGDATA_PATH))   # Multi-channel TIF
BINARY_DATA = BINARY_DATA .> 0 # Convert BINARY_DATA to Bool if needed (Ilastik often produces 0..255)

# Separate channels in OG_DATA if channels are interleaved in the 3rd dimension
# (Assumes shape is (height, width, colorSlices).
OG_DATA_CH1 = OG_DATA[:, :, 1:CHANNELS:end]
OG_DATA_CH2 = CHANNELS > 1 ? OG_DATA[:, :, 2:CHANNELS:end] : nothing
OG_DATA_CH3 = CHANNELS > 2 ? OG_DATA[:, :, 3:CHANNELS:end] : nothing
OG_DATA_CH4 = CHANNELS > 3 ? OG_DATA[:, :, 4:CHANNELS:end] : nothing

## Morphological Dilation & Fill Holes
kernel = ones(Bool, 3, 3, 3)
BINARY_DATA = ImageMorphology.dilate(BINARY_DATA, kernel)
function fill_3d_holes(vol)
    invVol = .!vol
    se = ImageMorphology.strel(ones(Bool, 3, 3, 3))
    labeled = ImageMorphology.label_components(invVol, se) 
    outside_label = labeled[1,1,1]
    filled = copy(vol)
    filled[labeled .!= outside_label] .= true
    return filled
end
BINARY_DATA = fill_3d_holes(BINARY_DATA)

## Keep Only Largest Connected Component
se = ImageMorphology.strel(ones(Bool, 3, 3, 3))
labeled = ImageMorphology.label_components(BINARY_DATA, se)
label_counts = Dict{Int, Int}()
for i in labeled
    if i != 0  # Skip background (0)
        label_counts[i] = get(label_counts, i, 0) + 1
    end
end
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
BINARY_DATA .= (labeled .== biggest_label)

## Multiply the Final Binary Mask by the Original Data
BINARY_DATA_U16 = convert.(UInt16, BINARY_DATA)
NEW_CH1 = BINARY_DATA_U16 .* OG_DATA_CH1
NEW_CH2 = CHANNELS > 1 ? BINARY_DATA_U16 .* OG_DATA_CH2 : nothing
NEW_CH3 = CHANNELS > 2 ? BINARY_DATA_U16 .* OG_DATA_CH3 : nothing
NEW_CH4 = CHANNELS > 3 ? BINARY_DATA_U16 .* OG_DATA_CH4 : nothing

## Volume, Surface Area, and Principal Axis
MAG       = 10.0
BIN       = 1.0
Z_STEP    = 6.25
um2Pixel  = 6.25   # microscope scale
pixel2len = um2Pixel * BIN / MAG
pixel2area   = pixel2len^2
voxel2volume = pixel2len^2 * Z_STEP

voxelcount = count(BINARY_DATA) #Volume
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

## Save the Final Mask and Channel Data
#TODO, get this to better interact with imageJ's UInt16 expectations
 PATH = "./"
 save(joinpath(PATH, "newseg.tif"), BINARY_DATA)
 save(joinpath(PATH, "newdata_CH1.tif"), NEW_CH1)
 save(joinpath(PATH, "newdata_CH2.tif"), NEW_CH2)
 save(joinpath(PATH, "newdata_CH3.tif"), NEW_CH3)
 save(joinpath(PATH, "newdata_CH4.tif"), NEW_CH4)

# println("Processing complete. Files saved to $PATH")

## View, Plot 3D
figure = Figure()
ax = Axis3(figure[1, 1], 
perspectiveness=0.5,
azimuth=2.19,
elevation=0.57)
voxels!(ax, UInt8.(BINARY_DATA), color = :skyblue, alpha = 0.5)
display(figure)