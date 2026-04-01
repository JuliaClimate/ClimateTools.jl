# Maps

ClimateTools returns YAXArray/Cube values that can be plotted with your preferred library.

Example with Plots.jl:

```julia
using Plots

arr = Array(annualmax(cube))
heatmap(arr[:, :, 1], title="Annual Maximum (first time step)")
```
