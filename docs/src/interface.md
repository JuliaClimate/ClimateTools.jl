# Interface

## Merging ClimGrid type

Sometimes, the timeseries are split among multiple files (e.g. climate models outputs). To obtain the complete timeseries, you can `merge` 2 `ClimGrid`. The method is based on the merging of 2 `AxisArrays` and is overloaded for the `ClimGrid` type.

```julia-repl
C = merge(C1::ClimGrid, C2::ClimGrid)
```

## Operators

Basic operators are overloaded on `ClimGrid`.


[`minimum`](@ref)
[`maximum`](@ref)
[`std`](@ref)
[`var`](@ref)
