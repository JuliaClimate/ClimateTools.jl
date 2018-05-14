# Interface

## Merging ClimGrid type

Sometimes, the timeseries are split among multiple files (e.g. climate models outputs). To obtain the complete timeseries, you can [`merge`](@ref) 2 `ClimGrid`. The method is based on the merging of 2 `AxisArrays` and is overloaded for the `ClimGrid` type.

```julia
C = merge(C1::ClimGrid, C2::ClimGrid)
```

To merge multiple ClimGrid form an array of files, [`load`](@ref) has a method that accepts an array of files to merge.

## Operators

Basic operators are overloaded on `ClimGrid`.

[`mean`](@ref)
[`minimum`](@ref)
[`maximum`](@ref)
[`std`](@ref)
[`var`](@ref)
