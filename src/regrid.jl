"""
    Regridder(source::YAXArray, dest::YAXArray; method="bilinear",
              lonname_source=:longitude, latname_source=:latitude,
              lonname_dest=:longitude, latname_dest=:latitude)

Build a reusable regridding object inspired by xESMF: compute interpolation
weights once and apply them many times.

Supported methods:
- `"bilinear"` (alias: `"linear"`)
- `"nearest_s2d"` (alias: `"nearest"`)

Regridders can be persisted with `save_regridder` and `load_regridder` for
lightweight reuse across sessions.
"""
abstract type AbstractRegridWeights end

const REGRIDDER_SERIALIZATION_VERSION = 1

struct BilinearWeights <: AbstractRegridWeights
    i0::Vector{Int}
    i1::Vector{Int}
    j0::Vector{Int}
    j1::Vector{Int}
    w00::Vector{Float64}
    w10::Vector{Float64}
    w01::Vector{Float64}
    w11::Vector{Float64}
    nx_out::Int
    ny_out::Int
end

struct NearestWeights <: AbstractRegridWeights
    src_linear_index::Vector{Int}
    nx_out::Int
    ny_out::Int
end

struct Regridder{W<:AbstractRegridWeights}
    method::Symbol
    source_lon_name::Symbol
    source_lat_name::Symbol
    dest_lon_name::Symbol
    dest_lat_name::Symbol
    source_lon_flip::Bool
    source_lat_flip::Bool
    source_lon::Vector{Float64}
    source_lat::Vector{Float64}
    dest_lon::Vector{Float64}
    dest_lat::Vector{Float64}
    weights::W
end

function _normalize_regrid_method(method::String)
    m = lowercase(method)
    if m in ("linear", "bilinear")
        return :bilinear
    elseif m in ("nearest", "nearest_s2d")
        return :nearest_s2d
    elseif m == "idw"
        return :idw
    end
    error("Unsupported regridding method: $(method)")
end

function _resolve_dim_name(cube::YAXArray, requested::Symbol, candidates::Tuple)
    dim_names = Tuple(name.(cube.axes))
    if requested in dim_names
        return requested
    end

    for cand in candidates
        if cand in dim_names
            return cand
        end
    end

    error("Could not find dimension $(requested) or any of $(candidates) in cube axes $(dim_names)")
end

function _prepare_source_coordinate(coord::AbstractVector, label::Symbol)
    out = Float64.(collect(coord))
    isempty(out) && error("Source coordinate $(label) is empty.")

    if length(out) == 1
        return out, false
    end

    any(Base.diff(out) .== 0.0) && error("Source coordinate $(label) contains duplicate values.")

    if issorted(out)
        return out, false
    elseif issorted(out; rev=true)
        return reverse(out), true
    end

    error("Source coordinate $(label) must be monotonic (ascending or descending).")
end

function _bracket_position(x::Float64, grid::Vector{Float64})
    n = length(grid)
    n == 1 && return 1, 1, 0.0

    if x <= grid[1]
        return 1, 2, 0.0
    elseif x >= grid[end]
        return n - 1, n, 1.0
    end

    i0 = searchsortedlast(grid, x)
    i1 = i0 + 1
    denom = grid[i1] - grid[i0]
    t = denom == 0.0 ? 0.0 : (x - grid[i0]) / denom

    return i0, i1, clamp(t, 0.0, 1.0)
end

function _build_bilinear_weights(source_lon::Vector{Float64}, source_lat::Vector{Float64}, dest_lon::Vector{Float64}, dest_lat::Vector{Float64})
    nx_out = length(dest_lon)
    ny_out = length(dest_lat)
    nout = nx_out * ny_out

    i0 = Vector{Int}(undef, nout)
    i1 = Vector{Int}(undef, nout)
    j0 = Vector{Int}(undef, nout)
    j1 = Vector{Int}(undef, nout)
    w00 = Vector{Float64}(undef, nout)
    w10 = Vector{Float64}(undef, nout)
    w01 = Vector{Float64}(undef, nout)
    w11 = Vector{Float64}(undef, nout)

    k = 1
    for j in 1:ny_out
        y = dest_lat[j]
        jy0, jy1, ty = _bracket_position(y, source_lat)

        for i in 1:nx_out
            x = dest_lon[i]
            ix0, ix1, tx = _bracket_position(x, source_lon)

            i0[k] = ix0
            i1[k] = ix1
            j0[k] = jy0
            j1[k] = jy1

            if ix0 == ix1 && jy0 == jy1
                w00[k] = 1.0
                w10[k] = 0.0
                w01[k] = 0.0
                w11[k] = 0.0
            elseif ix0 == ix1
                w00[k] = 1.0 - ty
                w10[k] = 0.0
                w01[k] = ty
                w11[k] = 0.0
            elseif jy0 == jy1
                w00[k] = 1.0 - tx
                w10[k] = tx
                w01[k] = 0.0
                w11[k] = 0.0
            else
                w00[k] = (1.0 - tx) * (1.0 - ty)
                w10[k] = tx * (1.0 - ty)
                w01[k] = (1.0 - tx) * ty
                w11[k] = tx * ty
            end

            k += 1
        end
    end

    return BilinearWeights(i0, i1, j0, j1, w00, w10, w01, w11, nx_out, ny_out)
end

function _nearest_lon_index(source_lon::Vector{Float64}, x::Float64)
    best_i = 1
    best_d = abs(mod(source_lon[1] - x + 180.0, 360.0) - 180.0)

    for i in 2:length(source_lon)
        d = abs(mod(source_lon[i] - x + 180.0, 360.0) - 180.0)
        if d < best_d
            best_d = d
            best_i = i
        end
    end

    return best_i
end

function _build_nearest_weights(source_lon::Vector{Float64}, source_lat::Vector{Float64}, dest_lon::Vector{Float64}, dest_lat::Vector{Float64})
    nx_out = length(dest_lon)
    ny_out = length(dest_lat)
    nout = nx_out * ny_out
    src_linear_index = Vector{Int}(undef, nout)

    k = 1
    for j in 1:ny_out
        y = dest_lat[j]
        j_src = argmin(abs.(source_lat .- y))

        for i in 1:nx_out
            x = dest_lon[i]
            i_src = _nearest_lon_index(source_lon, x)
            src_linear_index[k] = i_src + (j_src - 1) * length(source_lon)
            k += 1
        end
    end

    return NearestWeights(src_linear_index, nx_out, ny_out)
end

function Regridder(source::YAXArray, dest::YAXArray; method::String="bilinear", lonname_source=:longitude, latname_source=:latitude, lonname_dest=:longitude, latname_dest=:latitude)
    source_lon_name = _resolve_dim_name(source, lonname_source, (:longitude, :lon, :x, :rlon))
    source_lat_name = _resolve_dim_name(source, latname_source, (:latitude, :lat, :y, :rlat))
    dest_lon_name = _resolve_dim_name(dest, lonname_dest, (:longitude, :lon, :x, :rlon))
    dest_lat_name = _resolve_dim_name(dest, latname_dest, (:latitude, :lat, :y, :rlat))

    source_lon, source_lon_flip = _prepare_source_coordinate(lookup(source, source_lon_name), source_lon_name)
    source_lat, source_lat_flip = _prepare_source_coordinate(lookup(source, source_lat_name), source_lat_name)

    dest_lon = Float64.(collect(lookup(dest, dest_lon_name)))
    dest_lat = Float64.(collect(lookup(dest, dest_lat_name)))

    method_symbol = _normalize_regrid_method(method)
    if method_symbol == :bilinear
        weights = _build_bilinear_weights(source_lon, source_lat, dest_lon, dest_lat)
    elseif method_symbol == :nearest_s2d
        weights = _build_nearest_weights(source_lon, source_lat, dest_lon, dest_lat)
    else
        error("Method $(method) is not supported by Regridder. Use bilinear/linear or nearest/nearest_s2d.")
    end

    return Regridder(
        method_symbol,
        source_lon_name,
        source_lat_name,
        dest_lon_name,
        dest_lat_name,
        source_lon_flip,
        source_lat_flip,
        source_lon,
        source_lat,
        dest_lon,
        dest_lat,
        weights,
    )
end

function _weights_payload(weights::BilinearWeights)
    return (
        type=:bilinear,
        i0=weights.i0,
        i1=weights.i1,
        j0=weights.j0,
        j1=weights.j1,
        w00=weights.w00,
        w10=weights.w10,
        w01=weights.w01,
        w11=weights.w11,
        nx_out=weights.nx_out,
        ny_out=weights.ny_out,
    )
end

function _weights_payload(weights::NearestWeights)
    return (
        type=:nearest_s2d,
        src_linear_index=weights.src_linear_index,
        nx_out=weights.nx_out,
        ny_out=weights.ny_out,
    )
end

function _weights_from_payload(payload)
    if payload.type == :bilinear
        return BilinearWeights(
            payload.i0,
            payload.i1,
            payload.j0,
            payload.j1,
            payload.w00,
            payload.w10,
            payload.w01,
            payload.w11,
            payload.nx_out,
            payload.ny_out,
        )
    elseif payload.type == :nearest_s2d
        return NearestWeights(payload.src_linear_index, payload.nx_out, payload.ny_out)
    end

    error("Unsupported serialized weights type: $(payload.type)")
end

function _regridder_payload(regridder::Regridder)
    return (
        version=REGRIDDER_SERIALIZATION_VERSION,
        method=regridder.method,
        source_lon_name=regridder.source_lon_name,
        source_lat_name=regridder.source_lat_name,
        dest_lon_name=regridder.dest_lon_name,
        dest_lat_name=regridder.dest_lat_name,
        source_lon_flip=regridder.source_lon_flip,
        source_lat_flip=regridder.source_lat_flip,
        source_lon=regridder.source_lon,
        source_lat=regridder.source_lat,
        dest_lon=regridder.dest_lon,
        dest_lat=regridder.dest_lat,
        weights=_weights_payload(regridder.weights),
    )
end

function _regridder_from_payload(payload)
    hasproperty(payload, :version) || error("Serialized regridder is missing a format version.")
    payload.version == REGRIDDER_SERIALIZATION_VERSION || error("Unsupported serialized regridder version: $(payload.version)")

    return Regridder(
        payload.method,
        payload.source_lon_name,
        payload.source_lat_name,
        payload.dest_lon_name,
        payload.dest_lat_name,
        payload.source_lon_flip,
        payload.source_lat_flip,
        Float64.(collect(payload.source_lon)),
        Float64.(collect(payload.source_lat)),
        Float64.(collect(payload.dest_lon)),
        Float64.(collect(payload.dest_lat)),
        _weights_from_payload(payload.weights),
    )
end

function _same_grid(a::AbstractVector{Float64}, b::AbstractVector{Float64}; atol::Float64=1e-10)
    length(a) == length(b) || return false
    return all(isapprox.(a, b; atol=atol, rtol=0.0))
end

function _source_grid_state(cube::YAXArray, regridder::Regridder)
    lon_name = _resolve_dim_name(cube, regridder.source_lon_name, (:longitude, :lon, :x, :rlon))
    lat_name = _resolve_dim_name(cube, regridder.source_lat_name, (:latitude, :lat, :y, :rlat))

    source_lon, source_lon_flip = _prepare_source_coordinate(lookup(cube, lon_name), lon_name)
    source_lat, source_lat_flip = _prepare_source_coordinate(lookup(cube, lat_name), lat_name)

    _same_grid(source_lon, regridder.source_lon) || error("Source longitude coordinates do not match the Regridder source grid.")
    _same_grid(source_lat, regridder.source_lat) || error("Source latitude coordinates do not match the Regridder source grid.")

    return lon_name, lat_name, source_lon_flip, source_lat_flip
end

@inline function _is_missing_or_nan(x)
    if ismissing(x)
        return true
    end
    if x isa AbstractFloat
        return isnan(x)
    end
    return false
end

function _weighted_value(vals::NTuple{4,Any}, weights::NTuple{4,Float64}; skipna::Bool=false, na_thres::Float64=1.0)
    valid = ntuple(i -> !_is_missing_or_nan(vals[i]), 4)

    if !skipna
        all(valid) || return NaN
        return sum(Float64(vals[i]) * weights[i] for i in 1:4)
    end

    totalw = sum(weights)
    validw = sum(valid[i] ? weights[i] : 0.0 for i in 1:4)
    totalw == 0.0 && return NaN

    missing_ratio = 1.0 - validw / totalw
    if validw == 0.0 || missing_ratio > na_thres
        return NaN
    end

    return sum(valid[i] ? Float64(vals[i]) * weights[i] : 0.0 for i in 1:4) / validw
end

function _apply_regridder!(xout, xin; regridder::Regridder, source_lon_flip::Bool=false, source_lat_flip::Bool=false, skipna::Bool=false, na_thres::Float64=1.0)
    field = xin
    if source_lon_flip
        field = reverse(field, dims=1)
    end
    if source_lat_flip
        field = reverse(field, dims=2)
    end

    if regridder.weights isa BilinearWeights
        w = regridder.weights
        k = 1
        for j in 1:w.ny_out
            for i in 1:w.nx_out
                vals = (
                    field[w.i0[k], w.j0[k]],
                    field[w.i1[k], w.j0[k]],
                    field[w.i0[k], w.j1[k]],
                    field[w.i1[k], w.j1[k]],
                )
                ww = (w.w00[k], w.w10[k], w.w01[k], w.w11[k])
                xout[i, j] = _weighted_value(vals, ww; skipna=skipna, na_thres=na_thres)
                k += 1
            end
        end
    elseif regridder.weights isa NearestWeights
        w = regridder.weights
        k = 1
        for j in 1:w.ny_out
            for i in 1:w.nx_out
                v = field[w.src_linear_index[k]]
                xout[i, j] = _is_missing_or_nan(v) ? NaN : Float64(v)
                k += 1
            end
        end
    else
        error("Unsupported regridder weights type: $(typeof(regridder.weights))")
    end

    return nothing
end

"""
    regrid(cube::YAXArray, regridder::Regridder; skipna=false, na_thres=1.0)

Apply a precomputed regridder to a cube. Extra non-spatial dimensions are
preserved. The source lon/lat coordinates must match the grid used to build or
load the regridder.
"""
function regrid(cube::YAXArray, regridder::Regridder; skipna::Bool=false, na_thres::Float64=1.0)
    (0.0 <= na_thres <= 1.0) || error("na_thres must be in [0, 1].")

    source_lon_name, source_lat_name, source_lon_flip, source_lat_flip = _source_grid_state(cube, regridder)

    indims = InDims(string(source_lon_name), string(source_lat_name))
    outdims = OutDims(
        Dim{regridder.dest_lon_name}(regridder.dest_lon),
        Dim{regridder.dest_lat_name}(regridder.dest_lat),
    )

    return mapCube(
        _apply_regridder!,
        cube;
        regridder=regridder,
        source_lon_flip=source_lon_flip,
        source_lat_flip=source_lat_flip,
        skipna=skipna,
        na_thres=na_thres,
        indims=indims,
        outdims=outdims,
        nthreads=Threads.nthreads(),
    )
end

(regridder::Regridder)(cube::YAXArray; kwargs...) = regrid(cube, regridder; kwargs...)

"""
    save_regridder(path, regridder)

Persist a regridder and its precomputed weights to a lightweight binary file for
reuse across sessions in the same ClimateTools/Julia environment.
"""
function save_regridder(path::AbstractString, regridder::Regridder)
    open(path, "w") do io
        Serialization.serialize(io, _regridder_payload(regridder))
    end

    return path
end

"""
    load_regridder(path) -> Regridder

Load a regridder saved with `save_regridder`.
"""
function load_regridder(path::AbstractString)
    payload = open(path, "r") do io
        Serialization.deserialize(io)
    end

    return _regridder_from_payload(payload)
end

"""
    regrid_cube(source::YAXArray, dest::YAXArray; kwargs...) -> YAXArray

Compatibility wrapper around the Regridder workflow.
"""
function regrid_cube(source::YAXArray, dest::YAXArray; lonname_source=:longitude, latname_source=:latitude, lonname_dest=:longitude, latname_dest=:latitude, method::String="linear", skipna::Bool=false, na_thres::Float64=1.0)
    regridder = Regridder(
        source,
        dest;
        method=method,
        lonname_source=lonname_source,
        latname_source=latname_source,
        lonname_dest=lonname_dest,
        latname_dest=latname_dest,
    )

    return regrid(source, regridder; skipna=skipna, na_thres=na_thres)
end

"""
    idw_griddata(pts, vals, londest2, latdest2; k=8, p=2)

Simple inverse-distance-weighted (IDW) scattered interpolator.
"""
function idw_griddata(pts::AbstractMatrix, vals::AbstractVector, londest2, latdest2; k::Int=8, p::Real=2)
    valid = .!isnan.(vals)
    pts_valid = pts[valid, :]
    vals_valid = vals[valid]

    if isempty(vals_valid)
        return fill(NaN, size(londest2))
    end

    nsrc = size(pts_valid, 1)
    nx, ny = size(londest2)
    out = Array{Float64}(undef, nx, ny)

    lon_src = pts_valid[:, 1]
    lat_src = pts_valid[:, 2]

    use_kdtree = true

    if use_kdtree
        try
            pts_base = hcat(lon_src, lat_src)
            pts_minus = hcat(lon_src .- 360.0, lat_src)
            pts_plus = hcat(lon_src .+ 360.0, lat_src)
            pts_all = vcat(pts_base, pts_minus, pts_plus)
            kdt = NearestNeighbors.KDTree(transpose(pts_all))

            nbase = size(pts_base, 1)
            map_index(idx) = mod1(idx, nbase)

            for j in 1:ny
                for i in 1:nx
                    xo = londest2[i, j]
                    yo = latdest2[i, j]
                    ksel = min(k, nbase)
                    inds, dists = NearestNeighbors.knn(kdt, [xo, yo], ksel, true)
                    dsel = sqrt.(dists)

                    hit = findfirst(x -> x < 1e-12, dsel)
                    if hit !== nothing
                        out[i, j] = vals_valid[map_index(inds[hit])]
                        continue
                    end

                    orig_inds = map(map_index, inds)
                    vsel = vals_valid[orig_inds]
                    w = 1.0 ./ (dsel .^ p)
                    wsum = sum(w)
                    out[i, j] = wsum == 0.0 ? NaN : sum(w .* vsel) / wsum
                end
            end

            return out
        catch
            use_kdtree = false
        end
    end

    for j in 1:ny
        for i in 1:nx
            xo = londest2[i, j]
            yo = latdest2[i, j]

            dx = mod.(lon_src .- xo .+ 180.0, 360.0) .- 180.0
            dy = lat_src .- yo
            d = sqrt.(dx .^ 2 .+ dy .^ 2)

            hit = findfirst(x -> x < 1e-12, d)
            if hit !== nothing
                out[i, j] = vals_valid[hit]
                continue
            end

            ksel = min(k, nsrc)
            inds = partialsortperm(d, 1:ksel)
            dsel = d[inds]
            vsel = vals_valid[inds]

            w = 1.0 ./ (dsel .^ p)
            wsum = sum(w)
            out[i, j] = wsum == 0.0 ? NaN : sum(w .* vsel) / wsum
        end
    end

    return out
end

"""
    rotated_to_geographic(rlon, rlat, grid_north_longitude, grid_north_latitude)

Convert rotated-pole coordinates (rlon, rlat) to geographic longitude/latitude.
"""
function rotated_to_geographic(rlon::AbstractArray, rlat::AbstractArray, grid_north_longitude::Real, grid_north_latitude::Real)
    deg2rad = pi / 180.0
    phi_r = rlat .* deg2rad
    lambda_r = rlon .* deg2rad
    beta = grid_north_latitude * deg2rad
    lambda_p = grid_north_longitude * deg2rad

    sin_phi = sin(beta) .* sin(phi_r) .+ cos(beta) .* cos(phi_r) .* cos(lambda_r)
    phi = asin.(sin_phi)

    y = cos(phi_r) .* sin(lambda_r)
    x = cos(beta) .* sin(phi_r) .- sin(beta) .* cos(phi_r) .* cos(lambda_r)
    lambda = lambda_p .+ atan.(y, x)

    lat = phi ./ deg2rad
    lon = lambda ./ deg2rad
    lon = mod.(lon .+ 180.0, 360.0) .- 180.0

    return lon, lat
end

function _destination_grid(londest, latdest)
    if ndims(londest) == 1 && ndims(latdest) == 1
        return ndgrid(Float64.(collect(londest)), Float64.(collect(latdest)))
    elseif ndims(londest) == 2 && ndims(latdest) == 2
        size(londest) == size(latdest) || error("2D destination lon/lat arrays must have the same shape.")
        return Float64.(londest), Float64.(latdest)
    end

    error("Destination grid must be either two 1D vectors or two 2D arrays.")
end

function _is_separable_grid(lon2d::AbstractArray, lat2d::AbstractArray; atol::Float64=1e-8)
    m, n = size(lon2d)
    cols_equal = all(Base.maximum(abs.(lon2d[:, j] .- lon2d[:, 1])) < atol for j in 1:n)
    rows_equal = all(Base.maximum(abs.(lat2d[i, :] .- lat2d[1, :])) < atol for i in 1:m)
    return cols_equal && rows_equal
end

function _regrid_curvilinear_field(loncurv::AbstractArray, latcurv::AbstractArray, data::AbstractArray, londest, latdest; method::String="linear")
    size(loncurv) == size(latcurv) || error("loncurv and latcurv must have the same shape.")
    size(loncurv) == size(data) || error("Curvilinear coordinates and data must have the same shape.")

    method_symbol = _normalize_regrid_method(method)
    londest2, latdest2 = _destination_grid(londest, latdest)
    data_float = Float64.(coalesce.(data, NaN))

    if method_symbol == :bilinear && _is_separable_grid(loncurv, latcurv)
        xg1 = vec(Float64.(loncurv[:, 1]))
        yg1 = vec(Float64.(latcurv[1, :]))

        if length(xg1) > 1 && Base.diff(xg1)[1] < 0.0
            xg1 = reverse(xg1)
            data_float = reverse(data_float, dims=1)
        end
        if length(yg1) > 1 && Base.diff(yg1)[1] < 0.0
            yg1 = reverse(yg1)
            data_float = reverse(data_float, dims=2)
        end

        Interpolations.deduplicate_knots!([xg1, yg1], move_knots=true)
        itp = interpolate((xg1, yg1), data_float, Gridded(Linear()))
        etpf = extrapolate(itp, Interpolations.Flat())

        out = similar(londest2, Float64)
        for j in axes(out, 2)
            for i in axes(out, 1)
                out[i, j] = etpf(londest2[i, j], latdest2[i, j])
            end
        end
        return out
    end

    pts = hcat(vec(Float64.(loncurv)), vec(Float64.(latcurv)))
    vals = vec(data_float)

    if method_symbol == :nearest_s2d
        return idw_griddata(pts, vals, londest2, latdest2; k=1, p=2)
    end

    return idw_griddata(pts, vals, londest2, latdest2; k=8, p=2)
end

"""
    regrid_curvilinear_to_regular(loncurv, latcurv, data, londest, latdest; method="linear")

Regrid a field from a curvilinear lon/lat grid to a regular destination grid.
"""
function regrid_curvilinear_to_regular(loncurv::AbstractArray, latcurv::AbstractArray, data::AbstractArray, londest, latdest; method::String="linear")
    return _regrid_curvilinear_field(loncurv, latcurv, data, londest, latdest; method=method)
end

"""
    regrid_rotated_curvilinear_to_regular(lonrot, latrot, data, londest, latdest;
                                          grid_north_longitude,
                                          grid_north_latitude,
                                          method="linear")

Regrid a field from a rotated curvilinear grid to a regular destination grid.
"""
function regrid_rotated_curvilinear_to_regular(lonrot::AbstractArray, latrot::AbstractArray, data::AbstractArray, londest, latdest; grid_north_longitude, grid_north_latitude, method::String="linear")
    lon_geo, lat_geo = rotated_to_geographic(lonrot, latrot, grid_north_longitude, grid_north_latitude)
    return _regrid_curvilinear_field(lon_geo, lat_geo, data, londest, latdest; method=method)
end

"""
    regrid_curvilinear_to_regular(source::YAXArray, dest::YAXArray; kwargs...)

YAXArrays interface for curvilinear-to-regular regridding of 2D fields.
"""
function regrid_curvilinear_to_regular(source::YAXArray, dest::YAXArray; lonname_source=:longitude, latname_source=:latitude, lonname_dest=:longitude, latname_dest=:latitude, method::String="linear")
    source_lon_name = _resolve_dim_name(source, lonname_source, (:longitude, :lon, :x, :rlon))
    source_lat_name = _resolve_dim_name(source, latname_source, (:latitude, :lat, :y, :rlat))
    dest_lon_name = _resolve_dim_name(dest, lonname_dest, (:longitude, :lon, :x, :rlon))
    dest_lat_name = _resolve_dim_name(dest, latname_dest, (:latitude, :lat, :y, :rlat))

    data = Array(source)
    ndims(data) == 2 || error("regrid_curvilinear_to_regular(source::YAXArray, dest::YAXArray) currently supports 2D source fields only.")

    lon_lookup = lookup(source, source_lon_name)
    lat_lookup = lookup(source, source_lat_name)
    if ndims(lon_lookup) == 1 && ndims(lat_lookup) == 1
        loncurv, latcurv = ndgrid(Float64.(collect(lon_lookup)), Float64.(collect(lat_lookup)))
    else
        loncurv = Float64.(Array(lon_lookup))
        latcurv = Float64.(Array(lat_lookup))
        size(loncurv) == size(data) || error("Curvilinear source longitude coordinates must match source data shape.")
        size(latcurv) == size(data) || error("Curvilinear source latitude coordinates must match source data shape.")
    end

    londest = Float64.(collect(lookup(dest, dest_lon_name)))
    latdest = Float64.(collect(lookup(dest, dest_lat_name)))

    out = _regrid_curvilinear_field(loncurv, latcurv, data, londest, latdest; method=method)
    return YAXArray((Dim{dest_lon_name}(londest), Dim{dest_lat_name}(latdest)), out)
end
