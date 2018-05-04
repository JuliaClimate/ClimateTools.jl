# Typical workflow for a climate scenario

Say you need to create a climate scenario over a region defined by the following polygon.

```julia-repl
julia> poly_reg = [[NaN -65 -80 -80 -65 -65];[NaN 42 42 52 52 42]]
2×6 Array{Float64,2}:
 NaN  -65.0  -80.0  -80.0  -65.0  -65.0
 NaN   42.0   42.0   52.0   52.0   42.0
```

First step is to extract data over this region, with the `nc2julia` function. Data can be downloaded over the ESGF nodes (e.g. https://esgf-node.llnl.gov/projects/esgf-llnl/). In the following example, a custom in-house simulation done at Ouranos is used, but any compliant CF file can be used.

```julia-repl
C = nc2julia("crcm_ouranos_bby_195306.nc", "tasmax", poly=poly_reg)
```

First step of verification is to map the time-mean data with the `mapclimgrid` function to see if there is something wrong.

```julia-repl
mapclimgrid(C)
```

Which should return the following map.

<p align="center">
  <img src="https://user-images.githubusercontent.com/3630311/39647931-c0926668-4fae-11e8-9bd9-64588229c49c.png?raw=true" width="771" height="388" alt="Temperature example"/>
</p>
