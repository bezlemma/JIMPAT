# TIFF Tag Constants (ID, Type) - Defined as UInt16
const TIFF_TAG_IMAGE_WIDTH = UInt16(256)
const TIFF_TAG_IMAGE_LENGTH = UInt16(257)
const TIFF_TAG_BITS_PER_SAMPLE = UInt16(258)
const TIFF_TAG_COMPRESSION = UInt16(259)
const TIFF_TAG_PHOTOMETRIC_INTERP = UInt16(262)
const TIFF_TAG_IMAGE_DESCRIPTION = UInt16(270)
const TIFF_TAG_STRIP_OFFSETS = UInt16(273)
const TIFF_TAG_SAMPLES_PER_PIXEL = UInt16(277)
const TIFF_TAG_ROWS_PER_STRIP = UInt16(278)
const TIFF_TAG_STRIP_BYTE_COUNTS = UInt16(279)
const TIFF_TAG_X_RESOLUTION = UInt16(282)
const TIFF_TAG_Y_RESOLUTION = UInt16(283)
const TIFF_TAG_PLANAR_CONFIGURATION = UInt16(284)
const TIFF_TAG_RESOLUTION_UNIT = UInt16(296)

# TIFF Data Type Constants - Defined as UInt16
const TIFF_TYPE_BYTE = UInt16(1)      # 8-bit unsigned integer
const TIFF_TYPE_ASCII = UInt16(2)     # 8-bit bytes w/ last char NUL
const TIFF_TYPE_SHORT = UInt16(3)     # 16-bit unsigned integer
const TIFF_TYPE_LONG = UInt16(4)      # 32-bit unsigned integer
const TIFF_TYPE_RATIONAL = UInt16(5)  # Two LONGs: numerator, denominator

# --- Helper Functions (Internal to the Module) ---

# Helper function to write a 12-byte IFD Tag Entry
function write_ifd_entry(io::IO, tag_id::UInt16, type::UInt16, count::UInt32, value_or_offset::UInt32)
    write(io, tag_id)          # Tag ID (2 bytes)
    write(io, type)            # Data Type (2 bytes)
    write(io, count)           # Number of values (4 bytes)
    write(io, value_or_offset) # Value or Offset (4 bytes)
end

# Helper to write RATIONAL (2x LONG) data and return offset
function write_rational_data(io::IO, numerator::UInt32, denominator::UInt32)
    offset = position(io)
    write(io, numerator)
    write(io, denominator)
    return UInt32(offset)
end

"""
    SaveTiffPlz(filename::String, A::AbstractArray{T, 5}) where T<:Unsigned

Saves a 5D array `A` (order x,y,c,z,t) as a TIFF hyperstack readable by ImageJ/FIJI.

Handles `UInt8` and `UInt16` data types. 12-bit data should be passed in a `UInt16` array.
Writes uncompressed TIFF files. Reduces dimensions of size 1 implicitly in the loop structure
but includes them in the ImageJ metadata tag.

Args:
    filename (String): The path to save the TIFF file.
    A (AbstractArray{T, 5}): The 5D image array (x, y, c, z, t). T must be UInt8 or UInt16.

Limitations:
    - Assumes Little Endian system.
    - Only writes uncompressed data.
    - 12-bit data is saved as 16-bit.
    - Minimal metadata beyond ImageJ hyperstack structure.
    - Basic error handling.
"""
function SaveTiffPlz(filename::String, A::AbstractArray{T, 5}) where T<:Unsigned
    # --- 1. Validation ---
    if !(T <: Union{UInt8, UInt16})
        error("Unsupported data type: $T. Only UInt8 and UInt16 are supported.")
    end
    if isempty(A)
        error("Input array is empty.")
    end

    # --- 2. Get Dimensions & Basic Sizes ---
    Nx, Ny, Nc, Nz, Nt = size(A)
    num_planes = Nc * Nz * Nt
    if num_planes == 0
        error("Input array has zero planes (Nc*Nz*Nt = 0).")
    end


    bits_per_sample = sizeof(T) * 8
    if bits_per_sample != 8 && bits_per_sample != 16
        error("Calculated bits per sample ($bits_per_sample) is not 8 or 16.")
    end
    bytes_per_pixel = sizeof(T)
    # Ensure calculations use large enough types to avoid overflow before final casting
    bytes_per_plane_calc::Int64 = Int64(Nx) * Int64(Ny) * Int64(bytes_per_pixel) # Use Int64 for calculation
    if bytes_per_plane_calc > typemax(UInt32)
         @warn "Plane size ($bytes_per_plane_calc bytes) exceeds typemax(UInt32). TIFFs > 4GB might have issues."
         # Consider erroring out or implementing TIFF 64-bit extensions (BigTIFF - not done here)
    end
    bytes_per_plane = UInt32(bytes_per_plane_calc) # Store as UInt32 for TIFF tags


    # --- 3. Prepare ImageJ Metadata String ---
    ij_metadata = @sprintf("ImageJ=1.53k\nimages=%d\nchannels=%d\nslices=%d\nframes=%d\nhyperstack=true\nmode=grayscale\n",
                           num_planes, Nc, Nz, Nt)
    ij_metadata_bytes = vcat(Vector{UInt8}(ij_metadata), [0x00]) # Add Null terminator
    ij_metadata_len = length(ij_metadata_bytes) # This is Int

    # --- 4. Calculate Sizes & Offsets --- ## CORRECTED SECTION ##
    # Count the tags *actually* written in each IFD:
    # Standard tags: Width, Length, BPS, Compr, Photo, StripOff, SPP, RPS, StripBC, Planar, ResUnit, XRes, YRes
    num_standard_tags_written = 13
    num_tags_first_written = num_standard_tags_written + 1 # Add ImageDescription for the first IFD
    num_tags_other_written = num_standard_tags_written

    bytes_per_ifd_entry = 12
    # Size = NumEntriesField (2) + Entries (N*12) + NextIFDOffset (4)
    # Use the corrected tag counts here:
    size_ifd_first::Int64 = 2 + num_tags_first_written * bytes_per_ifd_entry + 4
    size_ifd_other::Int64 = 2 + num_tags_other_written * bytes_per_ifd_entry + 4
    # Need Int64 for num_planes in case it's huge, and use max(0,...) in case num_planes=1
    total_ifd_block_size::Int64 = size_ifd_first + max(0, Int64(num_planes) - 1) * size_ifd_other # Use Int64

    # Determine starting offsets using Int64 for calculations to avoid premature overflow
    offset_first_ifd_calc::Int64 = 8
    offset_rational_data_calc::Int64 = offset_first_ifd_calc + total_ifd_block_size # This is now correct
    offset_xres_data_calc::Int64 = offset_rational_data_calc
    offset_yres_data_calc::Int64 = offset_xres_data_calc + 8
    offset_ij_metadata_calc::Int64 = offset_yres_data_calc + 8
    offset_first_image_plane_calc::Int64 = offset_ij_metadata_calc + ij_metadata_len # This is now correct
    if offset_first_image_plane_calc % 2 != 0 # Optional word alignment
        offset_first_image_plane_calc += 1
    end

    # Check if offsets exceed UInt32 max *before* casting
    # The final address needed is start of last plane + size of last plane
    final_offset_needed::Int64 = offset_first_image_plane_calc + max(0, Int64(num_planes)-1)*Int64(bytes_per_plane_calc)
    if final_offset_needed > typemax(UInt32)
         @warn "Calculated file size potentially exceeds 4GB (UInt32 max offset: $final_offset_needed > $(typemax(UInt32))). Standard TIFF may not support this."
         # Consider erroring or implementing BigTIFF
    end

    # Cast final offsets needed for tags to UInt32
    offset_first_ifd = UInt32(offset_first_ifd_calc)
    offset_xres_data = UInt32(offset_xres_data_calc)
    offset_yres_data = UInt32(offset_yres_data_calc)
    offset_ij_metadata = UInt32(offset_ij_metadata_calc)
    offset_first_image_plane = UInt32(offset_first_image_plane_calc)
    ## --- END OF CORRECTED SECTION 4 --- ##

    # --- 5. Write TIFF File ---
    try
        open(filename, "w") do io
            # --- Write Header (8 Bytes) ---
            write(io, UInt8('I'), UInt8('I')) # Byte order: Little Endian
            write(io, UInt16(42))             # TIFF Version
            write(io, offset_first_ifd)       # Offset to the first IFD (already UInt32)

            # --- Write IFDs Sequentially ---
            current_ifd_offset_calc::Int64 = offset_first_ifd_calc # Use Int64 for tracking position
            rational_offsets_written = false
            # xres_data_offset and yres_data_offset are already UInt32

            plane_idx = 0
            for t in 1:Nt, z in 1:Nz, c in 1:Nc
                plane_idx += 1
                is_first_plane = (plane_idx == 1)

                # Use corrected tag counts for current IFD
                num_tags_current = is_first_plane ? num_tags_first_written : num_tags_other_written
                size_ifd_current = is_first_plane ? size_ifd_first : size_ifd_other

                # Calculate offset for the *next* IFD using Int64 arithmetic
                offset_next_ifd_calc::Int64 = (plane_idx == num_planes) ? 0 : (current_ifd_offset_calc + size_ifd_current)

                # Cast next offset to UInt32 only when writing
                offset_next_ifd = UInt32(offset_next_ifd_calc)
                 if offset_next_ifd_calc > typemax(UInt32) && offset_next_ifd_calc != 0 # Check potential overflow just before write
                     @warn "Next IFD offset ($offset_next_ifd_calc) exceeds UInt32 max."
                 end


                # Calculate offset for this plane's pixel data using UInt32 arithmetic (as base offset is UInt32)
                plane_offset_increment = UInt32(plane_idx - 1) * bytes_per_plane # bytes_per_plane is UInt32
                offset_image_data = offset_first_image_plane + plane_offset_increment # Result is UInt32

                # --- Write IFD Header ---
                # Seek takes Int usually, use calculated Int64 offset. Check if it fits UInt32 for safety.
                if current_ifd_offset_calc > typemax(UInt32)
                    error("Current IFD offset exceeds UInt32 max, cannot write standard TIFF.")
                end
                seek(io, UInt32(current_ifd_offset_calc))
                write(io, UInt16(num_tags_current)) # Number of entries in this IFD

                # --- Write Tag Entries (Ensuring 5th arg is UInt32) ---
                write_ifd_entry(io, TIFF_TAG_IMAGE_WIDTH, TIFF_TYPE_LONG, UInt32(1), UInt32(Nx))
                write_ifd_entry(io, TIFF_TAG_IMAGE_LENGTH, TIFF_TYPE_LONG, UInt32(1), UInt32(Ny))
                write_ifd_entry(io, TIFF_TAG_BITS_PER_SAMPLE, TIFF_TYPE_SHORT, UInt32(1), UInt32(bits_per_sample))
                write_ifd_entry(io, TIFF_TAG_COMPRESSION, TIFF_TYPE_SHORT, UInt32(1), UInt32(1))
                write_ifd_entry(io, TIFF_TAG_PHOTOMETRIC_INTERP, TIFF_TYPE_SHORT, UInt32(1), UInt32(1))
                write_ifd_entry(io, TIFF_TAG_STRIP_OFFSETS, TIFF_TYPE_LONG, UInt32(1), offset_image_data) # Already UInt32
                write_ifd_entry(io, TIFF_TAG_SAMPLES_PER_PIXEL, TIFF_TYPE_SHORT, UInt32(1), UInt32(1))
                write_ifd_entry(io, TIFF_TAG_ROWS_PER_STRIP, TIFF_TYPE_LONG, UInt32(1), UInt32(Ny))
                write_ifd_entry(io, TIFF_TAG_STRIP_BYTE_COUNTS, TIFF_TYPE_LONG, UInt32(1), bytes_per_plane) # Already UInt32
                write_ifd_entry(io, TIFF_TAG_PLANAR_CONFIGURATION, TIFF_TYPE_SHORT, UInt32(1), UInt32(1))
                write_ifd_entry(io, TIFF_TAG_RESOLUTION_UNIT, TIFF_TYPE_SHORT, UInt32(1), UInt32(1))

                # --- Resolution Tags ---
                if !rational_offsets_written
                    current_pos_ifd = position(io) # Remember position within IFD entries
                    seek(io, offset_xres_data) # Rational data offset (already UInt32)
                    # write_rational_data returns UInt32 offset, store temporarily if needed for verification
                    local written_xres_offset = write_rational_data(io, UInt32(1), UInt32(1)) # 1/1 resolution
                    local written_yres_offset = write_rational_data(io, UInt32(1), UInt32(1)) # 1/1 resolution
                    # Optional: Verify against precalculated offsets
                    # if written_xres_offset != offset_xres_data || written_yres_offset != offset_yres_data
                    #    error("Rational data offset mismatch during write.")
                    # end
                    seek(io, current_pos_ifd) # Go back to writing IFD entries
                    rational_offsets_written = true
                end
                # Write offset to rational data (already UInt32)
                write_ifd_entry(io, TIFF_TAG_X_RESOLUTION, TIFF_TYPE_RATIONAL, UInt32(1), offset_xres_data)
                write_ifd_entry(io, TIFF_TAG_Y_RESOLUTION, TIFF_TYPE_RATIONAL, UInt32(1), offset_yres_data)


                # --- ImageDescription Tag (FIRST IFD ONLY) ---
                if is_first_plane
                    # Write offset to ImageJ metadata (already UInt32)
                    write_ifd_entry(io, TIFF_TAG_IMAGE_DESCRIPTION, TIFF_TYPE_ASCII, UInt32(ij_metadata_len), offset_ij_metadata)
                end

                # --- Write Offset to Next IFD ---
                write(io, offset_next_ifd) # Already UInt32

                # --- Update Offset for Next Loop Iteration ---
                current_ifd_offset_calc += size_ifd_current # Use Int64 tracking variable

            end # End of IFD writing loop

            # --- Write Tag Data (ASCII Image Description) ---
            # Rational data was already written above.
            seek(io, offset_ij_metadata) # Already UInt32
            write(io, ij_metadata_bytes)

            # --- Write Pixel Data ---
            plane_idx = 0
            for t in 1:Nt, z in 1:Nz, c in 1:Nc
                plane_idx += 1
                # Recalculate offset exactly as before to ensure correctness
                plane_offset_increment = UInt32(plane_idx - 1) * bytes_per_plane
                offset_image_data = offset_first_image_plane + plane_offset_increment # Should be UInt32
                seek(io, offset_image_data)

                # Using a view avoids allocating memory for the plane slice
                plane_data = @view A[:, :, c, z, t]
                write(io, plane_data)
            end # End of pixel data writing loop

        end # Closes io implicitly
    catch e
        println(stderr, "Error writing TIFF file '$filename':")
        showerror(stderr, e)
        # If using Logging module, use @error instead
        println(stderr)
        # Clean up partially written file
        isfile(filename) && rm(filename; force=true)
        rethrow(e)
    finally
        # Optional: Add cleanup here if needed, though `open...do` handles file closing
    end

    println("Successfully saved hyperstack to $filename")
    return nothing
end