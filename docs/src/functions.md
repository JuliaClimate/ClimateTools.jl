# Index

```@docs
inpoly(p, poly::Matrix)
inpolygrid(lon::AbstractArray{N, 2} where N, lat::AbstractArray{N,2} where N, poly::AbstractArray{N,2} where N)
interp_climgrid(A::ClimGrid, B::ClimGrid; method::String="linear", min=[], max=[])
interp_climgrid(A::ClimGrid, lon::AbstractArray{N, 1} where N, lat::AbstractArray{N, 1} where N; method::String="linear", min=[], max=[])
meshgrid{T}(vx::AbstractVector{T}, vy::AbstractVector{T})
```
