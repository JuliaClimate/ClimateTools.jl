module ClimateToolsGeoMakieExt

using ClimateTools
using GeoMakie
using YAXArrays
using DimensionalData
using Dates
using Statistics

const Makie = GeoMakie.Makie
const DEFAULT_SOURCE_PROJ = "+proj=longlat +datum=WGS84"
const ENSEMBLE_DIM_CANDIDATES = (:member, :members, :ensemble, :scenario, :realization, :realisation, :run)

function _is_missing_or_nan(value)
    return ismissing(value) || (value isa Number && isnan(value))
end

function _finite_values(values)
    output = Float64[]
    for value in values
        if !_is_missing_or_nan(value)
            numeric = Float64(value)
            isfinite(numeric) && push!(output, numeric)
        end
    end
    return output
end

function _float_vector(values)
    return [ismissing(value) ? NaN : Float64(value) for value in values]
end

function _axis_names(cube::YAXArray)
    return Tuple(name.(cube.axes))
end

function _resolve_time_axis(cube::YAXArray)
    return ClimateTools._time_dim_symbol(cube)
end

function _resolve_selectors(cube::YAXArray, selectors::NamedTuple)
    indices = Any[Colon() for _ in 1:ndims(cube)]
    names = _axis_names(cube)

    for (key, value) in pairs(selectors)
        position = findfirst(==(key), names)
        position === nothing && error("Selector dimension $(key) not found in cube axes $(names).")
        indices[position] = value
    end

    return indices
end

function _apply_selectors(cube::YAXArray, selectors::NamedTuple)
    isempty(selectors) && return cube
    indices = _resolve_selectors(cube, selectors)
    return view(cube, indices...)
end

function _selectors_with(selectors::NamedTuple, key::Symbol, value)
    haskey(selectors, key) && error("Selector dimension $(key) cannot also be used as a varying plotting dimension.")
    return (; selectors..., key => value)
end

function _remaining_dim_names(cube::YAXArray, excluded::Tuple)
    return [dim for dim in _axis_names(cube) if !(dim in excluded)]
end

function _reduce_values(values, reducer)
    finite = _finite_values(values)
    isempty(finite) && return NaN
    if reducer === nothing
        length(finite) == 1 && return first(finite)
        error("Multiple values remain after selection. Provide a `reducer` or further selectors.")
    end
    return reducer(finite)
end

function _slice_nonspatial_dims(cube::YAXArray, lon_name::Symbol, lat_name::Symbol; dim=nothing, index=1)
    names = _axis_names(cube)
    lon_pos = findfirst(==(lon_name), names)
    lat_pos = findfirst(==(lat_name), names)
    extra_names = [names[i] for i in eachindex(names) if i != lon_pos && i != lat_pos]

    isempty(extra_names) && return cube

    slice_dim = if isnothing(dim)
        if length(extra_names) == 1
            first(extra_names)
        elseif :time in extra_names
            :time
        else
            error("geomap requires a single non-spatial dimension or an explicit `dim` keyword. Remaining dimensions: $(extra_names)")
        end
    else
        Symbol(dim)
    end

    position = findfirst(==(slice_dim), names)
    position === nothing && error("Slice dimension $(slice_dim) not found in cube axes $(names).")

    indices = Any[Colon() for _ in 1:ndims(cube)]
    indices[position] = index
    sliced = view(cube, indices...)

    remaining = setdiff(_axis_names(sliced), (lon_name, lat_name))
    isempty(remaining) || error("geomap still has non-spatial dimensions after slicing: $(remaining). Use `selectors` to reduce them first.")

    return sliced
end

function _regular_field(cube::YAXArray; selectors=(;), dim=nothing, index=1)
    selected = _apply_selectors(cube, selectors)
    lon_name = ClimateTools._find_spatial_dim(selected, (:longitude, :lon, :x, :rlon))
    lat_name = ClimateTools._find_spatial_dim(selected, (:latitude, :lat, :y, :rlat))
    lon_name === nothing && error("Could not infer longitude axis for plotting.")
    lat_name === nothing && error("Could not infer latitude axis for plotting.")

    grid_mapping = ClimateTools._detect_grid_mapping(selected)
    if grid_mapping !== nothing && (lon_name == :rlon || lat_name == :rlat || lowercase(grid_mapping) in ("rotated_pole", "rotated_latitude_longitude"))
        error("This cube uses grid_mapping=$(repr(grid_mapping)) and should be plotted through `geomap(dataset, :varname)` or `geomapfacet(dataset, :varname)` so ClimateTools can recover geographic coordinates.")
    end

    sliced = _slice_nonspatial_dims(selected, lon_name, lat_name; dim=dim, index=index)
    names = _axis_names(sliced)
    lon_pos = findfirst(==(lon_name), names)
    lat_pos = findfirst(==(lat_name), names)
    field = Array(sliced)
    field2d = permutedims(field, (lon_pos, lat_pos))

    lon = Float64.(collect(lookup(sliced, lon_name)))
    lat = Float64.(collect(lookup(sliced, lat_name)))

    return (; lon=lon, lat=lat, field=field2d)
end

function _dataset_field(ds::Dataset, varname::Symbol; selectors=(;), dim=nothing, index=1)
    cube = ds[varname]
    selected = _apply_selectors(cube, selectors)
    lon_name = ClimateTools._find_spatial_dim(selected, (:longitude, :lon, :x, :rlon))
    lat_name = ClimateTools._find_spatial_dim(selected, (:latitude, :lat, :y, :rlat))
    lon_name === nothing && error("Could not infer longitude axis for plotting.")
    lat_name === nothing && error("Could not infer latitude axis for plotting.")

    sliced = _slice_nonspatial_dims(selected, lon_name, lat_name; dim=dim, index=index)
    names = _axis_names(sliced)
    lon_pos = findfirst(==(lon_name), names)
    lat_pos = findfirst(==(lat_name), names)
    field = Array(sliced)
    field2d = permutedims(field, (lon_pos, lat_pos))

    grid_mapping = ClimateTools._detect_grid_mapping(sliced)
    if grid_mapping === nothing
        lon = Float64.(collect(lookup(sliced, lon_name)))
        lat = Float64.(collect(lookup(sliced, lat_name)))
        return (; lon=lon, lat=lat, field=field2d)
    end

    lon2d, lat2d = ClimateTools._extract_geographic_coords(ds, sliced, grid_mapping)
    return (; lon=Float64.(lon2d), lat=Float64.(lat2d), field=field2d)
end

function _normalize_regular_grid(lon::AbstractVector, lat::AbstractVector, field::AbstractMatrix; lon_0::Real=0.0)
    lon_out = collect(lon)
    lat_out = collect(lat)
    field_out = copy(field)

    if length(lon_out) > 1 && issorted(lon_out; rev=true)
        lon_out = reverse(lon_out)
        field_out = reverse(field_out; dims=1)
    end

    if length(lat_out) > 1 && issorted(lat_out; rev=true)
        lat_out = reverse(lat_out)
        field_out = reverse(field_out; dims=2)
    end

    lon_span = Base.maximum(lon_out) - Base.minimum(lon_out)
    if lon_span > 300
        shifted = mod.(lon_out .- lon_0 .+ 180.0, 360.0) .- 180.0 .+ lon_0
        perm = sortperm(shifted)
        lon_out = shifted[perm]
        field_out = field_out[perm, :]
    end

    return (; lon=lon_out, lat=lat_out, field=field_out)
end

function _normalized_map_payload(payload; lon_0::Real=0.0)
    if payload.lon isa AbstractVector && payload.lat isa AbstractVector
        return _normalize_regular_grid(payload.lon, payload.lat, payload.field; lon_0=lon_0)
    end

    return payload
end

function _axis_extent(values)
    finite = vec(values)[isfinite.(vec(values))]
    isempty(finite) && error("Could not infer projection extent from empty or non-finite coordinates.")
    return Base.minimum(finite), Base.maximum(finite)
end

function _default_dest(lon, lat; lon_0::Real=0.0)
    lon_min, lon_max = _axis_extent(lon)
    lat_min, lat_max = _axis_extent(lat)
    lon_span = lon_max - lon_min
    lat_span = lat_max - lat_min

    if lon_span > 300 && lat_span > 120
        return "+proj=eqearth +lon_0=$(lon_0)"
    end

    return DEFAULT_SOURCE_PROJ
end

function _default_map_title(title, dim, index)
    title !== nothing && return title
    isnothing(dim) && return "ClimateTools map"
    return "$(dim) = $(index)"
end

function _panel_title(panel_titles, facetdim::Symbol, value, panel_index::Int)
    if panel_titles === nothing
        return string(facetdim, " = ", value)
    elseif panel_titles isa Function
        return panel_titles(value, panel_index)
    else
        return panel_titles[panel_index]
    end
end

function _colorrange_from_payloads(payloads)
    mins = Float64[]
    maxs = Float64[]

    for payload in payloads
        values = _finite_values(vec(payload.field))
        isempty(values) && continue
        push!(mins, Base.minimum(values))
        push!(maxs, Base.maximum(values))
    end

    isempty(mins) && return nothing
    return (Base.minimum(mins), Base.maximum(maxs))
end

function _resolve_dest(payload, dest, lon_0)
    isnothing(dest) || return dest
    return _default_dest(payload.lon, payload.lat; lon_0=lon_0)
end

function _default_dest_from_limits(lonlims::Tuple, latlims::Tuple; lon_0::Real=0.0)
    lon_span = lonlims[2] - lonlims[1]
    lat_span = latlims[2] - latlims[1]

    if lon_span > 300 && lat_span > 120
        return "+proj=eqearth +lon_0=$(lon_0)"
    end

    return DEFAULT_SOURCE_PROJ
end

function _expand_limit_pair(lo::Real, hi::Real; fraction::Real=0.0)
    if lo == hi
        pad = max(abs(Float64(lo)) * 0.01, 1e-6)
        return (Float64(lo) - pad, Float64(hi) + pad)
    end

    span = Float64(hi) - Float64(lo)
    pad = span * Float64(fraction)
    return (Float64(lo) - pad, Float64(hi) + pad)
end

function _payload_limits(payload; lon_0::Real=0.0, padding_fraction::Tuple{<:Real,<:Real}=(0.0, 0.0))
    normalized = _normalized_map_payload(payload; lon_0=lon_0)
    lon_min, lon_max = _axis_extent(normalized.lon)
    lat_min, lat_max = _axis_extent(normalized.lat)
    lonlims = _expand_limit_pair(lon_min, lon_max; fraction=padding_fraction[1])
    latlims = _expand_limit_pair(lat_min, lat_max; fraction=padding_fraction[2])
    return (lonlims, latlims)
end

function _shared_payload_limits(payloads; lon_0::Real=0.0, padding_fraction::Tuple{<:Real,<:Real}=(0.0, 0.0))
    lon_mins = Float64[]
    lon_maxs = Float64[]
    lat_mins = Float64[]
    lat_maxs = Float64[]

    for payload in payloads
        normalized = _normalized_map_payload(payload; lon_0=lon_0)
        lon_min, lon_max = _axis_extent(normalized.lon)
        lat_min, lat_max = _axis_extent(normalized.lat)
        push!(lon_mins, lon_min)
        push!(lon_maxs, lon_max)
        push!(lat_mins, lat_min)
        push!(lat_maxs, lat_max)
    end

    isempty(lon_mins) && return nothing

    lonlims = _expand_limit_pair(Base.minimum(lon_mins), Base.maximum(lon_maxs); fraction=padding_fraction[1])
    latlims = _expand_limit_pair(Base.minimum(lat_mins), Base.maximum(lat_maxs); fraction=padding_fraction[2])
    return (lonlims, latlims)
end

function _facet_limits(payloads; limits=nothing, shared_spatial_limits::Bool=false, fit_limits::Bool=true, lon_0::Real=0.0, padding_fraction::Tuple{<:Real,<:Real}=(0.0, 0.0))
    limits !== nothing && return limits
    shared_spatial_limits || return nothing
    fit_limits || return nothing
    return _shared_payload_limits(payloads; lon_0=lon_0, padding_fraction=padding_fraction)
end

function _shared_projection_flag(shared_projection::Bool, shared_dest)
    isnothing(shared_dest) && return shared_projection
    return Bool(shared_dest)
end

function _shared_dest_value(payloads; dest, lon_0::Real=0.0)
    if isnothing(dest)
        shared_limits = _shared_payload_limits(payloads; lon_0=lon_0)
        shared_limits === nothing && return nothing
        return _default_dest_from_limits(shared_limits[1], shared_limits[2]; lon_0=lon_0)
    elseif dest isa Function
        error("When `shared_projection=true`, `dest` must be a single projection string or `nothing`, not a function.")
    elseif dest isa AbstractVector || dest isa Tuple
        error("When `shared_projection=true`, `dest` must be a single projection string or `nothing`, not a collection of per-panel projections.")
    end

    return dest
end

function _panel_destinations(payloads, labels; dest=nothing, shared_projection::Bool=true, lon_0::Real=0.0)
    if shared_projection
        shared_dest = _shared_dest_value(payloads; dest=dest, lon_0=lon_0)
        return fill(shared_dest, length(payloads))
    end

    if isnothing(dest)
        return [nothing for _ in payloads]
    elseif dest isa Function
        return [dest(payloads[index], labels[index], index) for index in eachindex(payloads)]
    elseif dest isa AbstractVector || dest isa Tuple
        length(dest) == length(payloads) || error("Per-panel `dest` must have the same length as the number of facet panels.")
        return collect(dest)
    end

    return fill(dest, length(payloads))
end

function _resolved_limits(payload; limits=nothing, fit_limits::Bool=true, lon_0::Real=0.0, padding_fraction::Tuple{<:Real,<:Real}=(0.0, 0.0))
    limits !== nothing && return limits
    fit_limits || return nothing
    return _payload_limits(payload; lon_0=lon_0, padding_fraction=padding_fraction)
end

function _apply_axis_visibility!(ax, frame::Bool)
    frame && return ax
    Makie.hidedecorations!(ax)
    Makie.hidespines!(ax)
    return ax
end

function _draw_coastlines!(ax; color=:black, linewidth=1.0)
    Makie.lines!(ax, GeoMakie.coastlines(); color=color, linewidth=linewidth)
end

function _render_map_axis!(slot, payload; title=nothing, source=DEFAULT_SOURCE_PROJ, dest=nothing, lon_0::Real=0.0, colormap=:viridis, colorrange=nothing, coastlines::Bool=true, coastline_color=:black, coastline_width=1.0, frame::Bool=true, limits=nothing, fit_limits::Bool=true, padding_fraction::Tuple{<:Real,<:Real}=(0.0, 0.0), axis_kwargs=(;), surface_kwargs=(;))
    normalized = _normalized_map_payload(payload; lon_0=lon_0)
    dest_proj = _resolve_dest(normalized, dest, lon_0)
    resolved_limits = _resolved_limits(normalized; limits=limits, fit_limits=fit_limits, lon_0=lon_0, padding_fraction=padding_fraction)
    axis_options = isnothing(resolved_limits) ? axis_kwargs : merge((; limits=resolved_limits), axis_kwargs)
    ax = GeoMakie.GeoAxis(slot; source=source, dest=dest_proj, title=title, axis_options...)
    plotobj = if isnothing(colorrange)
        Makie.surface!(ax, normalized.lon, normalized.lat, normalized.field;
            shading=Makie.NoShading,
            colormap=colormap,
            surface_kwargs...,
        )
    else
        Makie.surface!(ax, normalized.lon, normalized.lat, normalized.field;
            shading=Makie.NoShading,
            colormap=colormap,
            colorrange=colorrange,
            surface_kwargs...,
        )
    end

    coastlines && _draw_coastlines!(ax; color=coastline_color, linewidth=coastline_width)
    _apply_axis_visibility!(ax, frame)
    return ax, plotobj
end

function _add_colorbar!(slot, plotobj; label=nothing, vertical=true, colorbar_kwargs=(;))
    if isnothing(label)
        return Makie.Colorbar(slot, plotobj; vertical=vertical, colorbar_kwargs...)
    end

    return Makie.Colorbar(slot, plotobj; label=label, vertical=vertical, colorbar_kwargs...)
end

function _single_map_layout(fig, colorbar::Bool, colorbar_position::Symbol)
    if !colorbar
        return (; axis_slot=fig[1, 1], colorbar_slot=nothing, vertical=true)
    elseif colorbar_position == :right
        return (; axis_slot=fig[1, 1], colorbar_slot=fig[1, 2], vertical=true)
    elseif colorbar_position == :bottom
        return (; axis_slot=fig[1, 1], colorbar_slot=fig[2, 1], vertical=false)
    else
        error("Unsupported colorbar_position $(colorbar_position). Supported positions are :right and :bottom.")
    end
end

function _facet_layout_slot(fig, panel_index::Int, ncols::Int)
    row = cld(panel_index, ncols)
    col = mod1(panel_index, ncols)
    return row, col, fig[row, col]
end

function _facet_colorbar_slot(fig, nrows::Int, ncols::Int, colorbar_position::Symbol)
    if colorbar_position == :right
        return fig[1:nrows, ncols + 1], true
    elseif colorbar_position == :bottom
        return fig[nrows + 1, 1:ncols], false
    else
        error("Unsupported colorbar_position $(colorbar_position). Supported positions are :right and :bottom.")
    end
end

function _facet_indices(cube::YAXArray, facetdim::Symbol; selectors=(;), indices=nothing, maxpanels=nothing)
    haskey(selectors, facetdim) && error("facetdim $(facetdim) cannot also appear in selectors.")
    selected = _apply_selectors(cube, selectors)
    names = _axis_names(selected)
    facetdim in names || error("facetdim $(facetdim) not found in cube axes $(names).")
    values = collect(lookup(selected, facetdim))
    chosen = isnothing(indices) ? collect(eachindex(values)) : collect(indices)

    if !isnothing(maxpanels)
        chosen = chosen[1:min(length(chosen), maxpanels)]
    end

    labels = [values[i] for i in chosen]
    return chosen, labels
end

function _map_payload(cube::YAXArray; selectors=(;), dim=nothing, index=1)
    return _regular_field(cube; selectors=selectors, dim=dim, index=index)
end

function _map_payload(ds::Dataset, varname::Symbol; selectors=(;), dim=nothing, index=1)
    return _dataset_field(ds, varname; selectors=selectors, dim=dim, index=index)
end

function _facet_payloads(cube::YAXArray; facetdim::Symbol, selectors=(;), indices=nothing, maxpanels=nothing)
    facet_indices, labels = _facet_indices(cube, facetdim; selectors=selectors, indices=indices, maxpanels=maxpanels)
    payloads = [_map_payload(cube; selectors=_selectors_with(selectors, facetdim, idx)) for idx in facet_indices]
    return payloads, labels
end

function _facet_payloads(ds::Dataset, varname::Symbol; facetdim::Symbol, selectors=(;), indices=nothing, maxpanels=nothing)
    facet_indices, labels = _facet_indices(ds[varname], facetdim; selectors=selectors, indices=indices, maxpanels=maxpanels)
    payloads = [_map_payload(ds, varname; selectors=_selectors_with(selectors, facetdim, idx)) for idx in facet_indices]
    return payloads, labels
end

function _timeseries_groupdim(selected::YAXArray, time_name::Symbol, mode::Symbol, groupdim)
    remaining = _remaining_dim_names(selected, (time_name,))

    if !isnothing(groupdim)
        return Symbol(groupdim)
    elseif mode == :stats && :stats in remaining
        return :stats
    else
        for candidate in ENSEMBLE_DIM_CANDIDATES
            candidate in remaining && return candidate
        end
        return nothing
    end
end

function _grouped_series_from_cube(cube::YAXArray; selectors=(;), groupdim, reducer=nothing)
    selected = _apply_selectors(cube, selectors)
    time_name = _resolve_time_axis(selected)
    names = _axis_names(selected)
    time_pos = findfirst(==(time_name), names)
    group_pos = findfirst(==(groupdim), names)
    group_pos === nothing && error("groupdim $(groupdim) not found in cube axes $(names).")
    arr = Array(selected)

    perm = (time_pos, group_pos, (i for i in 1:ndims(arr) if i != time_pos && i != group_pos)...)
    arranged = permutedims(arr, perm)
    matrix = reshape(arranged, size(arranged, 1), size(arranged, 2), :)
    output = fill(NaN, size(matrix, 1), size(matrix, 2))

    for time_index in axes(matrix, 1), group_index in axes(matrix, 2)
        output[time_index, group_index] = _reduce_values(view(matrix, time_index, group_index, :), reducer)
    end

    labels = collect(lookup(selected, groupdim))
    return _time_values(selected, time_name), labels, output, time_name
end

function _stats_indices(labels)
    label_map = Dict(lowercase(string(label)) => index for (index, label) in enumerate(labels))
    haskey(label_map, "min") || error("timeseriesplot(mode=:stats) expects a `:stats` dimension containing `min`, `mean`, and `max`.")
    haskey(label_map, "mean") || error("timeseriesplot(mode=:stats) expects a `:stats` dimension containing `min`, `mean`, and `max`.")
    haskey(label_map, "max") || error("timeseriesplot(mode=:stats) expects a `:stats` dimension containing `min`, `mean`, and `max`.")
    return label_map["min"], label_map["mean"], label_map["max"]
end

function _timeseries_figure(; figure_size=(900, 400), title="Time series", xlabel="Time", ylabel="Value")
    fig = Makie.Figure(size=figure_size)
    ax = Makie.Axis(fig[1, 1]; title=title, xlabel=xlabel, ylabel=ylabel)
    return fig, ax
end

function _plot_series_collection!(ax, x_values, labels, matrix; alpha=0.45, linewidth=1.5, legend=true)
    for (group_index, label) in enumerate(labels)
        Makie.lines!(ax, x_values, matrix[:, group_index]; label=string(label), alpha=alpha, linewidth=linewidth)
    end

    legend && Makie.axislegend(ax)
    return ax
end

function _ribbon_bounds(matrix; spread=:minmax)
    lower = fill(NaN, size(matrix, 1))
    upper = fill(NaN, size(matrix, 1))
    center = fill(NaN, size(matrix, 1))

    for row in axes(matrix, 1)
        finite = _finite_values(view(matrix, row, :))
        isempty(finite) && continue

        center[row] = mean(finite)
        if spread == :minmax
            lower[row] = minimum(finite)
            upper[row] = maximum(finite)
        elseif spread isa Tuple && length(spread) == 2
            lower[row] = quantile(finite, spread[1])
            upper[row] = quantile(finite, spread[2])
        else
            error("Unsupported spread $(spread). Use :minmax or a quantile tuple such as (0.1, 0.9).")
        end
    end

    return lower, center, upper
end

function _grouped_values_from_cube(cube::YAXArray; selectors=(;), groupdim)
    selected = _apply_selectors(cube, selectors)
    names = _axis_names(selected)
    group_pos = findfirst(==(groupdim), names)
    group_pos === nothing && error("groupdim $(groupdim) not found in cube axes $(names).")
    values = Array(selected)
    labels = collect(lookup(selected, groupdim))
    grouped = Vector{Vector{Float64}}(undef, length(labels))

    for idx in eachindex(labels)
        indices = Any[Colon() for _ in 1:ndims(selected)]
        indices[group_pos] = idx
        grouped[idx] = _finite_values(vec(values[indices...]))
    end

    return labels, grouped
end

function ClimateTools.geomap(cube::YAXArray; selectors=(;), dim=nothing, index=1, title=nothing, source=DEFAULT_SOURCE_PROJ, dest=nothing, lon_0::Real=0.0, colormap=:viridis, colorrange=nothing, coastlines::Bool=true, coastline_color=:black, coastline_width=1.0, colorbar::Bool=false, colorbar_label=nothing, colorbar_position::Symbol=:right, frame::Bool=true, limits=nothing, fit_limits::Bool=true, limit_padding::Tuple{<:Real,<:Real}=(0.0, 0.0), figure_size=(900, 500), axis_kwargs=(;), surface_kwargs=(;), colorbar_kwargs=(;))
    payload = _map_payload(cube; selectors=selectors, dim=dim, index=index)
    fig = Makie.Figure(size=figure_size)
    layout = _single_map_layout(fig, colorbar, colorbar_position)
    _, plotobj = _render_map_axis!(layout.axis_slot, payload;
        title=_default_map_title(title, dim, index),
        source=source,
        dest=dest,
        lon_0=lon_0,
        colormap=colormap,
        colorrange=colorrange,
        coastlines=coastlines,
        coastline_color=coastline_color,
        coastline_width=coastline_width,
        frame=frame,
        limits=limits,
        fit_limits=fit_limits,
        padding_fraction=limit_padding,
        axis_kwargs=axis_kwargs,
        surface_kwargs=surface_kwargs,
    )

    colorbar && _add_colorbar!(layout.colorbar_slot, plotobj; label=colorbar_label, vertical=layout.vertical, colorbar_kwargs=colorbar_kwargs)
    return fig
end

function ClimateTools.geomap(ds::Dataset, varname::Symbol; selectors=(;), dim=nothing, index=1, title=nothing, source=DEFAULT_SOURCE_PROJ, dest=nothing, lon_0::Real=0.0, colormap=:viridis, colorrange=nothing, coastlines::Bool=true, coastline_color=:black, coastline_width=1.0, colorbar::Bool=false, colorbar_label=nothing, colorbar_position::Symbol=:right, frame::Bool=true, limits=nothing, fit_limits::Bool=true, limit_padding::Tuple{<:Real,<:Real}=(0.0, 0.0), figure_size=(900, 500), axis_kwargs=(;), surface_kwargs=(;), colorbar_kwargs=(;))
    payload = _map_payload(ds, varname; selectors=selectors, dim=dim, index=index)
    fig = Makie.Figure(size=figure_size)
    layout = _single_map_layout(fig, colorbar, colorbar_position)
    _, plotobj = _render_map_axis!(layout.axis_slot, payload;
        title=_default_map_title(title, dim, index),
        source=source,
        dest=dest,
        lon_0=lon_0,
        colormap=colormap,
        colorrange=colorrange,
        coastlines=coastlines,
        coastline_color=coastline_color,
        coastline_width=coastline_width,
        frame=frame,
        limits=limits,
        fit_limits=fit_limits,
        padding_fraction=limit_padding,
        axis_kwargs=axis_kwargs,
        surface_kwargs=surface_kwargs,
    )

    colorbar && _add_colorbar!(layout.colorbar_slot, plotobj; label=colorbar_label, vertical=layout.vertical, colorbar_kwargs=colorbar_kwargs)
    return fig
end

function ClimateTools.geomapfacet(cube::YAXArray; facetdim, selectors=(;), indices=nothing, maxpanels=nothing, ncols=3, panel_titles=nothing, title=nothing, source=DEFAULT_SOURCE_PROJ, dest=nothing, lon_0::Real=0.0, colormap=:viridis, colorrange=nothing, shared_colorrange::Bool=true, shared_spatial_limits::Bool=false, shared_projection::Bool=true, shared_dest=nothing, coastlines::Bool=true, coastline_color=:black, coastline_width=1.0, colorbar::Bool=true, colorbar_label=nothing, colorbar_position::Symbol=:right, frame::Bool=true, limits=nothing, fit_limits::Bool=true, limit_padding::Tuple{<:Real,<:Real}=(0.0, 0.0), figure_size=nothing, axis_kwargs=(;), surface_kwargs=(;), colorbar_kwargs=(;))
    payloads, labels = _facet_payloads(cube; facetdim=Symbol(facetdim), selectors=selectors, indices=indices, maxpanels=maxpanels)
    panel_count = length(payloads)
    ncols = max(1, min(ncols, panel_count))
    nrows = cld(panel_count, ncols)
    resolved_size = isnothing(figure_size) ? (380 * ncols + (colorbar && colorbar_position == :right ? 120 : 0), 280 * nrows + (colorbar && colorbar_position == :bottom ? 110 : 0)) : figure_size
    resolved_colorrange = isnothing(colorrange) && shared_colorrange ? _colorrange_from_payloads(payloads) : colorrange
    resolved_limits = _facet_limits(payloads; limits=limits, shared_spatial_limits=shared_spatial_limits, fit_limits=fit_limits, lon_0=lon_0, padding_fraction=limit_padding)
    resolved_shared_projection = _shared_projection_flag(shared_projection, shared_dest)
    resolved_dests = _panel_destinations(payloads, labels; dest=dest, shared_projection=resolved_shared_projection, lon_0=lon_0)

    fig = Makie.Figure(size=resolved_size)
    title === nothing || Makie.Label(fig[0, 1:ncols], title, fontsize=22)
    first_plot = nothing

    for panel_index in eachindex(payloads)
        _, _, slot = _facet_layout_slot(fig, panel_index, ncols)
        _, plotobj = _render_map_axis!(slot, payloads[panel_index];
            title=_panel_title(panel_titles, Symbol(facetdim), labels[panel_index], panel_index),
            source=source,
            dest=resolved_dests[panel_index],
            lon_0=lon_0,
            colormap=colormap,
            colorrange=resolved_colorrange,
            coastlines=coastlines,
            coastline_color=coastline_color,
            coastline_width=coastline_width,
            frame=frame,
            limits=resolved_limits,
            fit_limits=fit_limits && isnothing(resolved_limits),
            padding_fraction=limit_padding,
            axis_kwargs=axis_kwargs,
            surface_kwargs=surface_kwargs,
        )
        first_plot === nothing && (first_plot = plotobj)
    end

    if colorbar && first_plot !== nothing
        cslot, vertical = _facet_colorbar_slot(fig, nrows, ncols, colorbar_position)
        _add_colorbar!(cslot, first_plot; label=colorbar_label, vertical=vertical, colorbar_kwargs=colorbar_kwargs)
    end

    return fig
end

function ClimateTools.geomapfacet(ds::Dataset, varname::Symbol; facetdim, selectors=(;), indices=nothing, maxpanels=nothing, ncols=3, panel_titles=nothing, title=nothing, source=DEFAULT_SOURCE_PROJ, dest=nothing, lon_0::Real=0.0, colormap=:viridis, colorrange=nothing, shared_colorrange::Bool=true, shared_spatial_limits::Bool=false, shared_projection::Bool=true, shared_dest=nothing, coastlines::Bool=true, coastline_color=:black, coastline_width=1.0, colorbar::Bool=true, colorbar_label=nothing, colorbar_position::Symbol=:right, frame::Bool=true, limits=nothing, fit_limits::Bool=true, limit_padding::Tuple{<:Real,<:Real}=(0.0, 0.0), figure_size=nothing, axis_kwargs=(;), surface_kwargs=(;), colorbar_kwargs=(;))
    payloads, labels = _facet_payloads(ds, varname; facetdim=Symbol(facetdim), selectors=selectors, indices=indices, maxpanels=maxpanels)
    panel_count = length(payloads)
    ncols = max(1, min(ncols, panel_count))
    nrows = cld(panel_count, ncols)
    resolved_size = isnothing(figure_size) ? (380 * ncols + (colorbar && colorbar_position == :right ? 120 : 0), 280 * nrows + (colorbar && colorbar_position == :bottom ? 110 : 0)) : figure_size
    resolved_colorrange = isnothing(colorrange) && shared_colorrange ? _colorrange_from_payloads(payloads) : colorrange
    resolved_limits = _facet_limits(payloads; limits=limits, shared_spatial_limits=shared_spatial_limits, fit_limits=fit_limits, lon_0=lon_0, padding_fraction=limit_padding)
    resolved_shared_projection = _shared_projection_flag(shared_projection, shared_dest)
    resolved_dests = _panel_destinations(payloads, labels; dest=dest, shared_projection=resolved_shared_projection, lon_0=lon_0)

    fig = Makie.Figure(size=resolved_size)
    title === nothing || Makie.Label(fig[0, 1:ncols], title, fontsize=22)
    first_plot = nothing

    for panel_index in eachindex(payloads)
        _, _, slot = _facet_layout_slot(fig, panel_index, ncols)
        _, plotobj = _render_map_axis!(slot, payloads[panel_index];
            title=_panel_title(panel_titles, Symbol(facetdim), labels[panel_index], panel_index),
            source=source,
            dest=resolved_dests[panel_index],
            lon_0=lon_0,
            colormap=colormap,
            colorrange=resolved_colorrange,
            coastlines=coastlines,
            coastline_color=coastline_color,
            coastline_width=coastline_width,
            frame=frame,
            limits=resolved_limits,
            fit_limits=fit_limits && isnothing(resolved_limits),
            padding_fraction=limit_padding,
            axis_kwargs=axis_kwargs,
            surface_kwargs=surface_kwargs,
        )
        first_plot === nothing && (first_plot = plotobj)
    end

    if colorbar && first_plot !== nothing
        cslot, vertical = _facet_colorbar_slot(fig, nrows, ncols, colorbar_position)
        _add_colorbar!(cslot, first_plot; label=colorbar_label, vertical=vertical, colorbar_kwargs=colorbar_kwargs)
    end

    return fig
end

function _time_values(cube::YAXArray, time_name::Symbol)
    values = collect(lookup(cube, time_name))
    if eltype(values) <: Dates.TimeType || eltype(values) <: Number
        return values
    end
    return collect(1:length(values))
end

function _series_from_cube(cube::YAXArray; selectors=(;), reducer=nothing)
    selected = _apply_selectors(cube, selectors)
    time_name = _resolve_time_axis(selected)
    names = _axis_names(selected)
    time_pos = findfirst(==(time_name), names)
    arr = Array(selected)

    if ndims(arr) == 1
        return _time_values(selected, time_name), _float_vector(arr), time_name
    end

    perm = (time_pos, (i for i in 1:ndims(arr) if i != time_pos)...)
    arranged = permutedims(arr, perm)
    matrix = reshape(arranged, size(arranged, 1), :)

    series = if size(matrix, 2) == 1
        _float_vector(vec(matrix[:, 1]))
    elseif reducer === nothing
        error("timeseriesplot requires `selectors` that isolate one series, or a `reducer` to collapse non-time dimensions.")
    else
        output = fill(NaN, size(matrix, 1))
        for i in axes(matrix, 1)
            vals = _finite_values(view(matrix, i, :))
            if !isempty(vals)
                output[i] = reducer(vals)
            end
        end
        output
    end

    return _time_values(selected, time_name), series, time_name
end

function ClimateTools.timeseriesplot(series::AbstractVector; x=nothing, label=nothing, title="Time series", xlabel="Time", ylabel="Value", figure_size=(900, 400), color=:steelblue)
    x_values = isnothing(x) ? collect(1:length(series)) : collect(x)
    y_values = _float_vector(series)

    fig, ax = _timeseries_figure(; figure_size=figure_size, title=title, xlabel=xlabel, ylabel=ylabel)
    Makie.lines!(ax, x_values, y_values; label=label, color=color)
    label === nothing || Makie.axislegend(ax)
    return fig
end

function ClimateTools.timeseriesplot(series::NamedTuple; x=nothing, title="Time series", xlabel="Time", ylabel="Value", figure_size=(900, 400))
    fig, ax = _timeseries_figure(; figure_size=figure_size, title=title, xlabel=xlabel, ylabel=ylabel)

    for (label, values) in pairs(series)
        x_values = isnothing(x) ? collect(1:length(values)) : collect(x)
        Makie.lines!(ax, x_values, _float_vector(values); label=string(label))
    end

    Makie.axislegend(ax)
    return fig
end

function ClimateTools.timeseriesplot(cube::YAXArray; selectors=(;), reducer=nothing, title="Time series", xlabel="Time", ylabel="Value", figure_size=(900, 400), label=nothing, groupdim=nothing, mode::Symbol=:auto, spread=:minmax, color=:steelblue, alpha=0.4, linewidth=1.5, legend::Bool=true)
    selected = _apply_selectors(cube, selectors)
    time_name = _resolve_time_axis(selected)
    resolved_groupdim = _timeseries_groupdim(selected, time_name, mode, groupdim)
    resolved_mode = mode

    if resolved_mode == :auto
        resolved_mode = isnothing(resolved_groupdim) ? :single : (resolved_groupdim == :stats ? :stats : :lines)
    end

    if resolved_mode == :single
        x_values, y_values, inferred_time_name = _series_from_cube(cube; selectors=selectors, reducer=reducer)
        return ClimateTools.timeseriesplot(y_values;
            x=x_values,
            label=label,
            title=title,
            xlabel=String(inferred_time_name),
            ylabel=ylabel,
            figure_size=figure_size,
            color=color,
        )
    end

    isnothing(resolved_groupdim) && error("Could not infer an ensemble or stats dimension. Pass `groupdim=...` or use selectors/reducer to obtain a single series.")
    x_values, labels, matrix, inferred_time_name = _grouped_series_from_cube(cube; selectors=selectors, groupdim=resolved_groupdim, reducer=reducer)
    fig, ax = _timeseries_figure(; figure_size=figure_size, title=title, xlabel=String(inferred_time_name), ylabel=ylabel)

    if resolved_mode == :lines
        _plot_series_collection!(ax, x_values, labels, matrix; alpha=alpha, linewidth=linewidth, legend=legend)
    elseif resolved_mode == :mean_ribbon
        lower, center, upper = _ribbon_bounds(matrix; spread=spread)
        Makie.band!(ax, x_values, lower, upper; color=(color, alpha))
        Makie.lines!(ax, x_values, center; color=color, linewidth=max(2.0, linewidth), label=label === nothing ? "mean" : label)
        legend && Makie.axislegend(ax)
    elseif resolved_mode == :stats
        min_idx, mean_idx, max_idx = _stats_indices(labels)
        lower = matrix[:, min_idx]
        center = matrix[:, mean_idx]
        upper = matrix[:, max_idx]
        Makie.band!(ax, x_values, lower, upper; color=(color, alpha))
        Makie.lines!(ax, x_values, center; color=color, linewidth=max(2.0, linewidth), label=label === nothing ? "mean" : label)
        legend && Makie.axislegend(ax)
    else
        error("Unsupported timeseriesplot mode $(resolved_mode). Supported modes are :single, :lines, :mean_ribbon, and :stats.")
    end

    return fig
end

function _values_for_stats(data::AbstractVector)
    return _finite_values(data)
end

function _values_for_stats(cube::YAXArray)
    return _finite_values(vec(Array(cube)))
end

function ClimateTools.statsplot(data::Union{AbstractVector, YAXArray}; kind::Symbol=:hist, bins=30, title="Statistical plot", xlabel="Value", ylabel="Count", figure_size=(900, 400), color=:steelblue, groupdim=nothing)
    if data isa YAXArray && !isnothing(groupdim)
        labels, grouped = _grouped_values_from_cube(data; groupdim=Symbol(groupdim))
        named = NamedTuple{Tuple(Symbol.(string.(labels)))}(Tuple(grouped))
        return ClimateTools.statsplot(named; kind=kind, bins=bins, title=title, xlabel=xlabel, ylabel=ylabel, figure_size=figure_size)
    end

    values = _values_for_stats(data)
    fig = Makie.Figure(size=figure_size)
    ax = Makie.Axis(fig[1, 1]; title=title, xlabel=xlabel, ylabel=ylabel)

    if kind == :hist
        Makie.hist!(ax, values; bins=bins, color=color)
    elseif kind == :boxplot
        Makie.boxplot!(ax, fill(1, length(values)), values; color=color)
        ax.xticks = ([1], ["data"])
    else
        error("Unsupported statsplot kind $(kind) for a single series. Supported kinds: :hist and :boxplot")
    end

    return fig
end

function ClimateTools.statsplot(data::NamedTuple; kind::Symbol=:hist, bins=30, title="Statistical plot", xlabel="Value", ylabel="Count", figure_size=(900, 400))
    fig = Makie.Figure(size=figure_size)
    ax = Makie.Axis(fig[1, 1]; title=title, xlabel=xlabel, ylabel=ylabel)
    labels = collect(keys(data))

    if kind == :hist
        for (label, values) in pairs(data)
            Makie.hist!(ax, _values_for_stats(values); bins=bins, label=string(label), alpha=0.5)
        end
        Makie.axislegend(ax)
    elseif kind == :boxplot
        for (index, label) in enumerate(labels)
            values = _values_for_stats(data[label])
            Makie.boxplot!(ax, fill(index, length(values)), values)
        end
        ax.xticks = (collect(1:length(labels)), string.(labels))
    elseif kind == :scatter
        length(labels) == 2 || error("statsplot(kind=:scatter) requires exactly two named series.")
        first_values = _values_for_stats(data[labels[1]])
        second_values = _values_for_stats(data[labels[2]])
        count = min(length(first_values), length(second_values))
        Makie.scatter!(ax, first_values[1:count], second_values[1:count])
        ax.xlabel = string(labels[1])
        ax.ylabel = string(labels[2])
    else
        error("Unsupported statsplot kind $(kind). Supported kinds are :hist, :boxplot, and :scatter")
    end

    return fig
end

end
