#=
Includes:

dilate_3d

erode_3d

extract_largest_component

function count_neighbors(vol, x, y, z)
    - Count neighbors in a 3x3x3 neighborhood

function is_border(vol, x, y, z)
    - Check if a voxel is at the border

=#
function dilate_3d(volume, iterations=1)
    result = copy(volume)
    temp = falses(size(volume))
    
    # Neighbors in 6-connectivity
    directions = [
        (1, 0, 0), (-1, 0, 0),
        (0, 1, 0), (0, -1, 0),
        (0, 0, 1), (0, 0, -1)
    ]
    
    for _ in 1:iterations
        # Copy current result to temp
        temp .= result
        
        # Perform dilation
        for i in 1:size(volume, 1)
            for j in 1:size(volume, 2)
                for k in 1:size(volume, 3)
                    # If this voxel is already set, skip
                    if result[i, j, k]
                        continue
                    end
                    
                    # Check if any neighbors are set
                    for (dx, dy, dz) in directions
                        nx, ny, nz = i + dx, j + dy, k + dz
                        
                        if 1 <= nx <= size(volume, 1) && 
                           1 <= ny <= size(volume, 2) && 
                           1 <= nz <= size(volume, 3) &&
                           temp[nx, ny, nz]
                            result[i, j, k] = true
                            break
                        end
                    end
                end
            end
        end
    end
    
    return result
end

function erode_3d(volume, iterations=1)
    result = copy(volume)
    temp = falses(size(volume))
    
    # Neighbors in 6-connectivity
    directions = [
        (1, 0, 0), (-1, 0, 0),
        (0, 1, 0), (0, -1, 0),
        (0, 0, 1), (0, 0, -1)
    ]
    
    for _ in 1:iterations
        # Copy current result to temp
        temp .= result
        
        # Perform erosion
        for i in 1:size(volume, 1)
            for j in 1:size(volume, 2)
                for k in 1:size(volume, 3)
                    # Skip if not set
                    if !result[i, j, k]
                        continue
                    end
                    
                    # Check if any neighbors are not set (in original)
                    for (dx, dy, dz) in directions
                        nx, ny, nz = i + dx, j + dy, k + dz
                        
                        if !(1 <= nx <= size(volume, 1) && 
                             1 <= ny <= size(volume, 2) && 
                             1 <= nz <= size(volume, 3) &&
                             temp[nx, ny, nz])
                            result[i, j, k] = false
                            break
                        end
                    end
                end
            end
        end
    end
    
    return result
end


# Function to count neighbors in a 3x3x3 neighborhood
function count_neighbors(vol, x, y, z)
    count = 0
    neighbors26 = [(dx, dy, dz) for dx in -1:1 for dy in -1:1 for dz in -1:1 if !(dx == 0 && dy == 0 && dz == 0)]
    for (dx, dy, dz) in neighbors26
        nx, ny, nz = x + dx, y + dy, z + dz
        if 1 <= nx <= size(vol, 1) && 
            1 <= ny <= size(vol, 2) && 
            1 <= nz <= size(vol, 3) && 
            vol[nx, ny, nz]
            count += 1
        end
    end
    return count
end

# Function to check if a voxel is at the border
function is_border(vol, x, y, z)
    neighbors6 = [(1, 0, 0), (-1, 0, 0), # Define the 6-connectivity neighborhood
                    (0, 1, 0), (0, -1, 0),
                    (0, 0, 1), (0, 0, -1)]

    for (dx, dy, dz) in neighbors6
        nx, ny, nz = x + dx, y + dy, z + dz
        if 1 <= nx <= size(vol, 1) && 
            1 <= ny <= size(vol, 2) && 
            1 <= nz <= size(vol, 3) && 
            !vol[nx, ny, nz]
            return true
        end
    end
    return false
end
    


# Modify the end of your skeletonization section:
function extract_largest_component(binary_volume)
    # Create empty volume to store the result
    result = falses(size(binary_volume))
    
    # If volume is empty, return empty result
    if !any(binary_volume)
        return result
    end
    
    # Store component labels
    labels = zeros(Int, size(binary_volume))
    
    # Possible moves in 6 directions (6-connectivity)
    directions = [
        (1, 0, 0), (-1, 0, 0),
        (0, 1, 0), (0, -1, 0),
        (0, 0, 1), (0, 0, -1)
    ]
    
    # Perform connected component labeling
    current_label = 0
    component_sizes = Int[]
    
    # Find all connected components
    for i in 1:lastindex(binary_volume, 1)
        for j in 1:lastindex(binary_volume, 2)
            for k in 1:lastindex(binary_volume, 3)
                # Skip if not a voxel or already labeled
                if !binary_volume[i, j, k] || labels[i, j, k] != 0
                    continue
                end
                
                # New component found
                current_label += 1
                push!(component_sizes, 0)
                
                # Start BFS from this voxel
                queue = [(i, j, k)]
                labels[i, j, k] = current_label
                component_sizes[current_label] += 1
                
                while !isempty(queue)
                    x, y, z = popfirst!(queue)
                    
                    # Check all neighbors
                    for (dx, dy, dz) in directions
                        nx, ny, nz = x + dx, y + dy, z + dz
                        
                        # Check if within bounds
                        if 1 <= nx <= size(binary_volume, 1) && 
                           1 <= ny <= size(binary_volume, 2) && 
                           1 <= nz <= size(binary_volume, 3)
                            
                            # If part of object and not yet labeled
                            if binary_volume[nx, ny, nz] && labels[nx, ny, nz] == 0
                                labels[nx, ny, nz] = current_label
                                component_sizes[current_label] += 1
                                push!(queue, (nx, ny, nz))
                            end
                        end
                    end
                end
            end
        end
    end
    
    # Find the largest component
    if isempty(component_sizes)
        return result  # No components found
    end
    
    largest_label = argmax(component_sizes)
    largest_size = component_sizes[largest_label]
    
    println("Found $(length(component_sizes)) connected components")
    println("Largest component has $largest_size voxels ($(round(100*largest_size/sum(component_sizes), digits=1))% of all voxels)")
    
    # Extract the largest component
    for i in 1:size(binary_volume, 1)
        for j in 1:size(binary_volume, 2)
            for k in 1:size(binary_volume, 3)
                if labels[i, j, k] == largest_label
                    result[i, j, k] = true
                end
            end
        end
    end
    
    return result
end


function nearest_index_naive(skeleton_pts, query)
    bestidx = 1
    bestdist = Inf
    @inbounds for i in 1:lastindex(skeleton_pts,1)
        dx = skeleton_pts[i,1] - query[1]
        dy = skeleton_pts[i,2] - query[2]
        dz = skeleton_pts[i,3] - query[3]
        dist2 = dx*dx + dy*dy + dz*dz
        if dist2 < bestdist
            bestdist = dist2
            bestidx  = i
        end
    end
    return bestidx
end