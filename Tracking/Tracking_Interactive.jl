##
using TiffImages
using ImageMorphology
using ImageFiltering
using Statistics
using DataFrames
using GLMakie
using Colors
using Base.Threads
using LsqFit
using Observables  # for creating and updating Observables

include("./Tracking_Funcs.jl")

## Load image stack
fpath = "C:/Users/uComp/Desktop/20250212/Test2.tif"
@info "Loading $fpath..."
imgStack = TiffImages.load(fpath) .|> Float16
imgStack = (imgStack .- minimum(imgStack)) ./ (maximum(imgStack) - minimum(imgStack))

## 1) Find bright circles in all frames -> returns a DataFrame of centroids
df = find_colloids_parallel(imgStack)
@info "Found $(nrow(df)) total objects."

df = link_trajectories(df)

# Compute an effective radius from the area (Area is in pixels => Area = π r^2)
df[!, :Radius] = sqrt.(df[!, :Area] ./ π)

maxFrame = maximum(df[!, :Frame])
tracks = groupby(df, :TrackID)

df_avg = combine(groupby(df, :Frame), :Radius => mean => :AvgRadius)

# --- Fit a power law (R = a * t^b) to average radius ---
model(t, p) = p[1] .* t .^ p[2]
p0 = [1.0, 1.0]
fit_result = curve_fit(model, df_avg.Frame, df_avg.AvgRadius, p0)
p_fit = fit_result.param
@info "Fitted parameters: a=$(p_fit[1]), b=$(p_fit[2])"

ts = range(minimum(df_avg.Frame), maximum(df_avg.Frame), length=100)
fitted_curve = model(ts, p_fit)

## ------------------- INTERACTIVE PLOTTING -------------------
fig = Figure()

# Left panel: image axis
ax_img = Axis(
    fig[1, 1],
    title = "Image / Track Annotation",
    aspect = DataAspect()
)

# Right panel: radius vs. time
ax_rad = Axis(
    fig[1, 2],
    title = "Radius over Time",
    xlabel = "Frame", ylabel = "Radius (pixels)",
    xscale = log10, yscale = log10
)

# Slider below the image axis
slider = Slider(fig[2, 1], range=1:maxFrame, width=Relative(1.0))
frame_index = slider.value

# Create an Observable image slice that depends on frame_index
img_obs = lift(frame_index) do f
    transpose(imgStack[:, :, f])
end

# Display the current image frame in the left axis
#img_obs = img_obs[end:-1:1, :];
image!(ax_img, img_obs, interpolate=false, colorrange=(0,1))

# Plot all tracks (gray) on the right axis
for trk in tracks
    trk_sorted = sort(trk, :Frame)
    lines!(ax_rad, trk_sorted.Frame, trk_sorted.Radius, color=(0.5,0.5,0.5,0.5))
end

# Plot the average radius (black) and the fitted curve (red dashed)
lines!(ax_rad, df_avg.Frame, df_avg.AvgRadius, color=:black, linewidth=4, label="Avg Radius")
lines!(ax_rad, ts, fitted_curve, color=:red, linestyle=:dash, linewidth=2, label="Power Law Fit")
axislegend(ax_rad)

#
# --- OBSERVABLES FOR HIGHLIGHT PLOTS ---
#
hl_x_track = Observable([NaN])
hl_y_track = Observable([NaN])
lines!(ax_img, hl_x_track, hl_y_track, color=:yellow, linewidth=2)

hl_x_pts = Observable([NaN])
hl_y_pts = Observable([NaN])
scatter!(ax_img, hl_x_pts, hl_y_pts, color=:yellow, markersize=8)

hl_x_rad = Observable([NaN])
hl_y_rad = Observable([NaN])
lines!(ax_rad, hl_x_rad, hl_y_rad, color=:blue, linewidth=3)

#
# --- FIND NEAREST CENTROID FUNCTION ---
#
function find_nearest_centroid(df::DataFrame, frame::Int, xq::Float64, yq::Float64)
    # Filter to only the rows from that frame
    mask = df[!, :Frame] .== frame
    sub = df[mask, :]

    # If no particles in this frame, return nothing
    if nrow(sub) == 0
        return nothing
    end

    # We'll compute distances to each centroid
    rows_sub = collect(eachrow(sub))
    dist = map(r -> ( (r[:CentroidX] - xq)^2 + (r[:CentroidY] - yq)^2) , rows_sub)
    idx_min = argmin(dist)

    # The row with minimal distance
    chosen_row = rows_sub[idx_min]
    chosen_track_id = chosen_row[:TrackID]
    return chosen_track_id, chosen_row
end

#
# --- CLICK HANDLER ON THE IMAGE AXIS ---
#
on(events(ax_img).mousebutton) do ev
    # We'll only proceed if this is a left-press inside the axis
    if ev.button == Mouse.left && ev.action == Mouse.press
        if is_mouseinside(ax_img)
            #
            # NEW: Clear old data so these Observables can accept new shapes
            #
            hl_x_track[] = [NaN]
            hl_y_track[] = [NaN]
            hl_x_pts[]   = [NaN]
            hl_y_pts[]   = [NaN]
            hl_x_rad[]   = [NaN]
            hl_y_rad[]   = [NaN]

            # Convert the clicked pixel coords to axis data coords
            xclick, yclick = mouseposition(ax_img)

            # Current frame from slider
            f = frame_index[]

            # Find the nearest centroid for that frame
            out = find_nearest_centroid(df, f, xclick, yclick)
            if out === nothing
                @warn "No nearest centroid found!!!"
                return
            end
            track_id, row = out

            # row is a DataFrameRow, so we can do row[:Column]
            cx = row[:CentroidX]
            cy = row[:CentroidY]
            @info "Selected TrackID=$track_id at x=$cx, y=$cy, frame=$f"

            # Now filter entire track
            track_mask = df[!, :TrackID] .== track_id
            mytrack = df[track_mask, :]
            mytrack = sort(mytrack, :Frame)

            # Update highlights in the image axis
            hl_x_track[] = mytrack[!, :CentroidX]
            hl_y_track[] = mytrack[!, :CentroidY]

            hl_x_pts[]   = mytrack[!, :CentroidX]
            hl_y_pts[]   = mytrack[!, :CentroidY]

            # Update highlights in the radius vs. time axis
            hl_x_rad[]   = mytrack[!, :Frame]
            hl_y_rad[]   = mytrack[!, :Radius]
        end
    end
end

# Finally, display the figure
display(fig)
