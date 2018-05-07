# Index

```@meta
CurrentModule = ClimateTools
```

```@autodocs
Modules = [ClimateTools]
Order   = [:function, :type]
```

<!-- ```@docs
inpoly(p, poly::Matrix)

inpolygrid(lon::AbstractArray{N, 2} where N, lat::AbstractArray{N,2} where N, poly::AbstractArray{N,2} where N)

interp_climgrid(A::ClimGrid, B::ClimGrid; method::String="linear", min=[], max=[])

interp_climgrid(A::ClimGrid, lon::AbstractArray{N, 1} where N, lat::AbstractArray{N, 1} where N; method::String="linear", min=[], max=[])

meshgrid{T}(vx::AbstractVector{T}, vy::AbstractVector{T})

nc2julia(file::String, variable::String; poly, start_date::Date, end_date::Date, data_units::String)

qqmap(obs::ClimGrid, ref::ClimGrid, fut::ClimGrid; method::String, detrend::Bool, window::Int, rankn::Int, thresnan::Float64, keep_original::Bool, interp = Linear(), extrap = Flat())

qqmap(obsvec::Array{N, 1} where N, refvec::Array{N, 1} where N, futvec::Array{N, 1} where N, datevec_obs, datevec_ref, datevec_fut; method::String, detrend::Bool, window::Int64, rankn::Int64, thresnan::Float64, keep_original::Bool, interp = Linear(), extrap = Flat())

shapefile_coords(poly::Shapefile.Polygon)

shapefile_coords_poly(poly::Shapefile.Polygon)

``` -->
