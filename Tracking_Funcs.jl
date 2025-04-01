## 
using TiffImages
using ImageMorphology, ImageFiltering
using Statistics, DataFrames
using GLMakie
using Colors
using Base.Threads
using LsqFit

# Helper function.
function ind2sub(a, i)
    i2s = CartesianIndices(a)
    i2s[i]
end

# center-of-mass refinement.
function refine_centroid(img, initial_position, r_win)
    y0, x0 = initial_position
    y_min = max(Int(floor(y0 - r_win)), 1)
    y_max = min(Int(ceil(y0 + r_win)), size(img, 1))
    x_min = max(Int(floor(x0 - r_win)), 1)
    x_max = min(Int(ceil(x0 + r_win)), size(img, 2))
    window = img[y_min:y_max, x_min:x_max]
    ys = collect(y_min:y_max)
    xs = collect(x_min:x_max)
    total_intensity = sum(window)
    if total_intensity == 0
        return (y0, x0)
    end
    weighted_y = 0.0
    weighted_x = 0.0
    for i in 1:lastindex(ys)
        for j in 1:lastindex(xs)
            weighted_y += window[i, j] * ys[i]
            weighted_x += window[i, j] * xs[j]
        end
    end
    refined_y = weighted_y / total_intensity
    refined_x = weighted_x / total_intensity
    return (refined_y, refined_x)
end

function is_gaussian_enough(img, center, r_win; sigma_guess=r_win/2, threshold=0.8)
    y0, x0 = center
    y_min = max(Int(floor(y0 - r_win)), 1)
    y_max = min(Int(ceil(y0 + r_win)), size(img, 1))
    x_min = max(Int(floor(x0 - r_win)), 1)
    x_max = min(Int(ceil(x0 + r_win)), size(img, 2))
    window = img[y_min:y_max, x_min:x_max]
    
    ys = collect(y_min:y_max)
    xs = collect(x_min:x_max)
    
    # Estimate amplitude A as the difference between the maximum and the median in the window.
    A = maximum(window) - median(window)
    B = median(window)
    
    # Construct the ideal Gaussian model.
    model = [A * exp(-((y - y0)^2 + (x - x0)^2) / (2 * sigma_guess^2)) + B for y in ys, x in xs]
    
    # Flatten the arrays for correlation computation.
    actual = vec(window)
    model_vec = vec(model)
    
    r = cor(actual, model_vec)
    return (r >= threshold), r
end

# Parallelized version of the colloid finder.
function find_colloids_parallel(imgStack; thresholdFactor=1.0, r_win=20, gaussianThreshold=0.9)
    nframes = size(imgStack, 3)
    # Preallocate an array to hold the results for each frame.
    results_by_frame = Vector{DataFrame}(undef, nframes)
    
    @threads for f in 1:nframes
        local_results = DataFrame(Frame=Int[], ObjectID=Int[],
                                  Area=Int[], CentroidX=Float64[], CentroidY=Float64[])
        frame = imgStack[:, :, f]
        frame_f32 = Float32.(frame)  # Convert Float16 to Float32 for processing.
        # Smooth the image with a 1-pixel Gaussian.
        smoothed = imfilter(frame_f32, Kernel.gaussian(1))
        
        t = mean(smoothed) + thresholdFactor * std(smoothed)
        binary = smoothed .> t
        
        # Identify connected components.
        labeled = label_components(binary)
        
        for lbl in unique(labeled)
            if lbl == 0  # Skip background.
                continue
            end
            inds = findall(labeled .== lbl)
            area = length(inds)
            sum_intensity = 0.0
            sum_y = 0.0
            sum_x = 0.0
            # Compute an intensity-weighted centroid over the connected region.
            for idx in inds
                i, j = Tuple(ind2sub(size(smoothed), idx))
                I_val = smoothed[i, j]
                sum_intensity += I_val
                sum_y += I_val * i
                sum_x += I_val * j
            end
            if sum_intensity == 0
                centroid_y = 0.0
                centroid_x = 0.0
            else
                centroid_y = sum_y / sum_intensity
                centroid_x = sum_x / sum_intensity
            end
            
            # Refine the centroid using a local center-of-mass refinement.
            refined_y, refined_x = refine_centroid(smoothed, (centroid_y, centroid_x), r_win)
            
            # Check if the candidate's intensity profile is well approximated by a Gaussian.
            good, corr = is_gaussian_enough(smoothed, (refined_y, refined_x), r_win;
                                              sigma_guess=r_win/2, threshold=gaussianThreshold)
            # Only add the candidate if the Gaussian match is good enough.
            if good
                push!(local_results, (Frame=f, ObjectID=0, Area=area,
                                        CentroidX=refined_x, CentroidY=refined_y))
            end
        end
        results_by_frame[f] = local_results
    end
    
    # Combine all local results.
    results = vcat(results_by_frame...)
    # Assign unique ObjectIDs sequentially.
    results[!, :ObjectID] = 1:nrow(results)
    return results
end

function link_trajectories(df; L=200.0, areaWeight=0.1)
    # Ensure the DataFrame is sorted by Frame.
    df = sort(df, :Frame)
    n = nrow(df)
    df[!, :TrackID] = Union{Missing, Int}[missing for _ in 1:nrow(df)]

    nextTrackID = 1
    # Process frame 1: assign each detection a new track.
    frame1_idx = findall(df.Frame .== 1)
    for idx in frame1_idx
        df[idx, :TrackID] = nextTrackID
        nextTrackID += 1
    end

    # Build dictionary for previous frame tracks:
    # key = TrackID, value = (CentroidX, CentroidY, Area)
    prev_tracks = Dict{Int, Tuple{Float64,Float64,Float64}}()
    for idx in frame1_idx
        track = df[idx, :TrackID]
        prev_tracks[track] = (df[idx, :CentroidX], df[idx, :CentroidY], df[idx, :Area])
    end

    maxFrame = maximum(df.Frame)
    for f in 2:maxFrame
        current_idx = findall(df.Frame .== f)
        # Build list of candidate links: (track, current_index, cost)
        candidate_links = []
        for (track, (x_prev, y_prev, area_prev)) in prev_tracks
            for idx in current_idx
                x_cur = df[idx, :CentroidX]
                y_cur = df[idx, :CentroidY]
                area_cur = df[idx, :Area]
                d = sqrt((x_cur - x_prev)^2 + (y_cur - y_prev)^2)
                if d > L
                    continue
                end
                cost = d + areaWeight * abs(area_cur - area_prev)
                push!(candidate_links, (track, idx, cost))
            end
        end
        # Greedy assignment: sort candidates by cost and assign if neither track nor detection has been used.
        sort!(candidate_links, by = x -> x[3])
        assigned_tracks = Set{Int}()
        assigned_current = Set{Int}()
        new_prev_tracks = Dict{Int, Tuple{Float64,Float64,Float64}}()
        for (track, idx, cost) in candidate_links
            if (track in assigned_tracks) || (idx in assigned_current)
                continue
            end
            # Assign detection idx to this track.
            df[idx, :TrackID] = track
            push!(assigned_tracks, track)
            push!(assigned_current, idx)
            new_prev_tracks[track] = (df[idx, :CentroidX], df[idx, :CentroidY], df[idx, :Area])
        end
        # For current detections not assigned, assign new tracks.
        for idx in current_idx
            if df[idx, :TrackID] === missing
                df[idx, :TrackID] = nextTrackID
                new_prev_tracks[nextTrackID] = (df[idx, :CentroidX], df[idx, :CentroidY], df[idx, :Area])
                nextTrackID += 1
            end
        end
        prev_tracks = new_prev_tracks
    end
    return df
end