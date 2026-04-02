function _normalize_xmap_dim(dim::Symbol)
    return dim
end

function _normalize_xmap_dim(dim::AbstractString)
    return Symbol(dim)
end

function _normalize_xmap_dims(dim)
    return (_normalize_xmap_dim(dim),)
end

function _normalize_xmap_dims(dims::Tuple)
    return Tuple(_normalize_xmap_dim(dim) for dim in dims)
end

function _normalize_xmap_dims(dims::AbstractVector)
    return Tuple(_normalize_xmap_dim(dim) for dim in dims)
end

function _normalize_xmap_output_axes()
    return ()
end

function _normalize_xmap_output_axes(output_axis)
    return (output_axis,)
end

function _normalize_xmap_output_axes(output_axes::Tuple)
    return output_axes
end

function _normalize_xmap_output_axes(output_axes::AbstractVector)
    return Tuple(output_axes)
end

function _xmap_window(cube::YAXArray, dims::Tuple{Vararg{Symbol}})
    return length(dims) == 1 ? (cube ⊘ first(dims)) : (cube ⊘ dims)
end

function _xmap_output_spec(output_axes::Tuple; outtype=nothing)
    if isnothing(outtype)
        return isempty(output_axes) ? XOutput() : XOutput(output_axes...)
    end

    return isempty(output_axes) ? XOutput(; outtype=outtype) : XOutput(output_axes...; outtype=outtype)
end

function _drop_reduced_xmap_dims(result::YAXArray, reduced_dims::Tuple{Vararg{Symbol}})
    drop_positions = Int[]

    for (index, axis) in pairs(result.axes)
        if Symbol(name(axis)) in reduced_dims && size(result, index) == 1
            push!(drop_positions, index)
        end
    end

    return isempty(drop_positions) ? result : dropdims(result; dims=Tuple(drop_positions))
end

function _xmap_call(kernel::Function, cube::YAXArray; reduced_dims, output_axes=(), inplace::Bool=true, function_args=(), function_kwargs=(;), outtype=nothing)
    return _xmap_call(kernel, (cube,); reduced_dims=reduced_dims, output_axes=output_axes, inplace=inplace, function_args=function_args, function_kwargs=function_kwargs, outtype=outtype)
end

function _xmap_call(kernel::Function, cubes::Tuple{Vararg{YAXArray}}; reduced_dims, output_axes=(), inplace::Bool=true, function_args=(), function_kwargs=(;), outtype=nothing)
    normalized_dims = _normalize_xmap_dims(reduced_dims)
    normalized_output_axes = _normalize_xmap_output_axes(output_axes)
    windowed_cubes = tuple((_xmap_window(cube, normalized_dims) for cube in cubes)...)

    result = xmap(
        kernel,
        windowed_cubes...;
        output=_xmap_output_spec(normalized_output_axes; outtype=outtype),
        inplace=inplace,
        function_args=function_args,
        function_kwargs=function_kwargs,
    )

    return _drop_reduced_xmap_dims(result, normalized_dims)
end