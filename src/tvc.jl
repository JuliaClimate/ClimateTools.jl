const DEFAULT_TVC_SCALES = [365, 183, 92, 46, 23, 12, 6, 3, 2]

struct TVCModel
    scales::Vector{Int}
    init_index::Int
    obs_mean::Vector{Float64}
    train_mean::Vector{Float64}
    transform::Matrix{Float64}
end

function _tvc_invalid(value)
    return ismissing(value) || (value isa Number && isnan(value))
end

function _tvc_float_vector(values)
    return Float64[_tvc_invalid(value) ? NaN : Float64(value) for value in values]
end

function _validate_tvc_scales(scales)
    isempty(scales) && error("scales must not be empty.")

    normalized = Int[scale for scale in scales]
    all(scale > 0 for scale in normalized) || error("scales must contain positive integers.")
    issorted(normalized; rev=true) || error("scales must be sorted in descending order.")

    return normalized
end

function _tvc_init_index(scales::AbstractVector{<:Integer})
    return sum(scales) - length(scales)
end

function _tvc_default_output(values, init_index::Int; keep_original::Bool=false)
    output = fill(NaN, length(values))

    if keep_original && init_index < length(values)
        raw_values = _tvc_float_vector(values)
        output[(init_index + 1):end] .= raw_values[(init_index + 1):end]
    end

    return output
end

function _tvc_rolling_mean(values, window::Int)
    n = length(values)
    prefix_sum = zeros(Float64, n + 1)
    prefix_count = zeros(Int, n + 1)

    for i in eachindex(values)
        prefix_sum[i + 1] = prefix_sum[i]
        prefix_count[i + 1] = prefix_count[i]

        value = values[i]
        if isfinite(value)
            prefix_sum[i + 1] += value
            prefix_count[i + 1] += 1
        end
    end

    smoothed = fill(NaN, n)
    for i in 1:n
        first_index = max(1, i - window + 1)
        total = prefix_sum[i + 1] - prefix_sum[first_index]
        count = prefix_count[i + 1] - prefix_count[first_index]
        if count > 0
            smoothed[i] = total / count
        end
    end

    return smoothed
end

function _tvc_decompose(values, scales)
    data = _tvc_float_vector(values)
    n = length(data)
    level_count = length(scales)

    filtered = fill(NaN, n, level_count)
    residual = copy(data)

    for (level, scale) in pairs(scales)
        smoothed = _tvc_rolling_mean(residual, scale)
        next_residual = fill(NaN, n)

        for i in eachindex(residual)
            if isfinite(residual[i]) && isfinite(smoothed[i])
                next_residual[i] = residual[i] - smoothed[i]
            end
        end

        filtered[:, level] .= smoothed
        filtered[isnan.(next_residual), level] .= NaN
        residual = next_residual
    end

    return hcat(filtered, residual)
end

function _tvc_truncate(matrix::AbstractMatrix, init_index::Int)
    if init_index >= size(matrix, 1)
        return matrix[1:0, :]
    end

    return matrix[(init_index + 1):end, :]
end

function _tvc_valid_row_mask(matrix::AbstractMatrix)
    mask = BitVector(undef, size(matrix, 1))

    for row in axes(matrix, 1)
        mask[row] = all(isfinite, @view matrix[row, :])
    end

    return mask
end

function _tvc_stats(matrix::AbstractMatrix)
    mask = _tvc_valid_row_mask(matrix)
    count(mask) >= 2 || return nothing

    filtered = matrix[mask, :]
    mean_vector = vec(mean(filtered, dims=1))
    centered = filtered .- permutedims(mean_vector)
    covariance = centered' * centered / (size(filtered, 1) - 1)

    return mean_vector, covariance
end

function _tvc_matrix_power(covariance::AbstractMatrix, exponent::Real; eig_floor::Real)
    eig_floor > 0 || error("eig_floor must be positive.")

    symmetric_covariance = Hermitian((covariance + covariance') / 2)
    decomposition = eigen(symmetric_covariance)
    scale = max(Base.maximum(abs.(decomposition.values)), 1.0)
    regularized = max.(decomposition.values, eig_floor * scale)

    return decomposition.vectors * Diagonal(regularized .^ exponent) * decomposition.vectors'
end

function _tvc_apply_components(model::TVCModel, matrix::AbstractMatrix)
    shift = model.obs_mean .- model.train_mean
    shifted = matrix .+ permutedims(shift)
    shifted_mean = vec(mean(shifted, dims=1))
    centered = shifted .- permutedims(shifted_mean)

    return permutedims(shifted_mean) .+ centered * model.transform
end

"""
    fit_tvc(obs::AbstractVector, raw_train::AbstractVector; scales=DEFAULT_TVC_SCALES, eig_floor::Real=1e-8)

Fit a Time Variability Correction model following Shao et al. (2024),
DOI 10.1029/2023MS003640.

The fitted model stores the multiscale mean and covariance transform estimated
from observations and the raw training series. Use [`apply_tvc`] to correct a
validation or projection series.
"""
function fit_tvc(obs::AbstractVector, raw_train::AbstractVector; scales=DEFAULT_TVC_SCALES, eig_floor::Real=1e-8)
    length(obs) == length(raw_train) || error("obs and raw_train must have the same length.")

    normalized_scales = _validate_tvc_scales(scales)
    init_index = _tvc_init_index(normalized_scales)
    init_index < length(obs) || error("Series length must exceed the TVC warm-up period of $(init_index) samples.")

    obs_matrix = _tvc_truncate(_tvc_decompose(obs, normalized_scales), init_index)
    train_matrix = _tvc_truncate(_tvc_decompose(raw_train, normalized_scales), init_index)

    obs_stats = _tvc_stats(obs_matrix)
    train_stats = _tvc_stats(train_matrix)
    isnothing(obs_stats) && error("obs does not contain enough valid samples after TVC decomposition.")
    isnothing(train_stats) && error("raw_train does not contain enough valid samples after TVC decomposition.")

    obs_mean, obs_covariance = obs_stats
    train_mean, train_covariance = train_stats
    transform = _tvc_matrix_power(train_covariance, -0.5; eig_floor=eig_floor) * _tvc_matrix_power(obs_covariance, 0.5; eig_floor=eig_floor)

    return TVCModel(normalized_scales, init_index, obs_mean, train_mean, transform)
end

"""
    apply_tvc(model::TVCModel, raw_val::AbstractVector; keep_original::Bool=false)

Apply a fitted Time Variability Correction model to a validation or projection
time series.

The returned series keeps the same length as `raw_val`. The leading warm-up
segment required by the multiscale rolling decomposition is filled with `NaN`.
If `keep_original=true`, rows after the warm-up segment that cannot be corrected
because of invalid values are copied from the raw series.
"""
function apply_tvc(model::TVCModel, raw_val::AbstractVector; keep_original::Bool=false)
    output = _tvc_default_output(raw_val, model.init_index; keep_original=keep_original)
    model.init_index < length(raw_val) || return output

    truncated = _tvc_truncate(_tvc_decompose(raw_val, model.scales), model.init_index)
    mask = _tvc_valid_row_mask(truncated)
    any(mask) || return output

    corrected_components = _tvc_apply_components(model, truncated[mask, :])
    output_tail = @view output[(model.init_index + 1):end]
    output_tail[mask] .= vec(sum(corrected_components, dims=2))

    return output
end

"""
    tvc(obs::AbstractVector, raw_train::AbstractVector, raw_val::AbstractVector; scales=DEFAULT_TVC_SCALES, eig_floor::Real=1e-8, keep_original::Bool=false)

Fit and apply Time Variability Correction in one step.
"""
function tvc(obs::AbstractVector, raw_train::AbstractVector, raw_val::AbstractVector; scales=DEFAULT_TVC_SCALES, eig_floor::Real=1e-8, keep_original::Bool=false)
    model = fit_tvc(obs, raw_train; scales=scales, eig_floor=eig_floor)
    return apply_tvc(model, raw_val; keep_original=keep_original)
end

function _tvc_kernel(dataout, obsvec, raw_trainvec, raw_valvec; scales, eig_floor::Real, keep_original::Bool)
    init_index = _tvc_init_index(scales)

    corrected = try
        tvc(obsvec, raw_trainvec, raw_valvec; scales=scales, eig_floor=eig_floor, keep_original=keep_original)
    catch
        _tvc_default_output(raw_valvec, init_index; keep_original=keep_original)
    end

    dataout .= corrected
end

"""
    tvc(obs::YAXArray, raw_train::YAXArray, raw_val::YAXArray; scales=DEFAULT_TVC_SCALES, eig_floor::Real=1e-8, keep_original::Bool=false, drop_leapday::Bool=true)

Apply Time Variability Correction of Shao et al. (2024), DOI
10.1029/2023MS003640, to every grid-cell time series of a YAXArray cube.

This method corrects covariance across and between multiple time scales while
preserving the time-event sequence of the validation series. The correction is
applied with `mapCube` over the time dimension. The output retains the
validation time axis, with the leading TVC warm-up segment filled with `NaN`.
"""
function tvc(obs::YAXArray, raw_train::YAXArray, raw_val::YAXArray; scales=DEFAULT_TVC_SCALES, eig_floor::Real=1e-8, keep_original::Bool=false, drop_leapday::Bool=true)
    normalized_scales = _validate_tvc_scales(scales)
    obs_timedim = _time_dim_symbol(obs)
    train_timedim = _time_dim_symbol(raw_train)
    val_timedim = _time_dim_symbol(raw_val)

    if drop_leapday
        obs = drop29thfeb(obs, dimt=obs_timedim)
        raw_train = drop29thfeb(raw_train, dimt=train_timedim)
        raw_val = drop29thfeb(raw_val, dimt=val_timedim)
    end

    length(lookup(obs, obs_timedim)) == length(lookup(raw_train, train_timedim)) || error("obs and raw_train must have the same time-axis length after preprocessing.")

    obs_indims = InDims(string(obs_timedim))
    train_indims = InDims(string(train_timedim))
    val_indims = InDims(string(val_timedim))
    outdims = OutDims(Dim{:time}(lookup(raw_val, val_timedim)))

    return mapCube(
        _tvc_kernel,
        (obs, raw_train, raw_val),
        indims=(obs_indims, train_indims, val_indims),
        outdims=outdims,
        scales=normalized_scales,
        eig_floor=eig_floor,
        keep_original=keep_original,
        nthreads=Threads.nthreads(),
    )
end