# Interface

## Merging ClimGrid type

Sometimes, the timeseries are split among multiple files (e.g. climate models outputs). To obtain the complete timeseries, you can [`merge`](@ref) 2 `ClimGrid`. The method is based on `AxisArrays` merging method and is overloaded for the `ClimGrid` type.

```julia
C = merge(C1::ClimGrid, C2::ClimGrid)
```

To merge multiple ClimGrid form an array of files, [`load`](@ref) has a method that accepts an array of files to merge.

## Operators

Basic statistical functions are overloaded on `ClimGrid`.

[`mean`](@ref)
[`minimum`](@ref)
[`maximum`](@ref)
[`std`](@ref)
[`var`](@ref)

Basic arithmetic operators are also loaded.

```julia
D = C + 2.0 # will add 2.0 to all elements of C
D = C::ClimGrid - A::ClimGrid # subtract A from C (useful for climatological difference between a future and historical period
D = C / A # Ratio of 2 ClimGrids
```

## Exporting

Exporting a `ClimGrid` to disk to a netCDF format can be done with the `write` function.

```julia
write(C::ClimGrid, filename::String)
```
