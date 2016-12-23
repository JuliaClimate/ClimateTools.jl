
<a id='ClimateTools.jl-documentation-1'></a>

# ClimateTools.jl documentation

- [ClimateTools.jl documentation](index.md#ClimateTools.jl-documentation-1)
    - [Installation](index.md#Installation-1)
    - [Functions](index.md#Functions-1)
    - [Index](index.md#Index-1)


<a id='Installation-1'></a>

## Installation


This package is registered in `METADATA.jl` and so can be installed using `Pkg.add`


```julia
Pkg.add("ClimateTools")
```


For the latest version, checkout master branch


```julia
Pkg.checkout("ClimateTools")
```


<a id='Functions-1'></a>

## Functions

<a id='ClimateTools.frostdays-Tuple{Array{Float64,1},StepRange{Date,Base.Dates.Day}}' href='#ClimateTools.frostdays-Tuple{Array{Float64,1},StepRange{Date,Base.Dates.Day}}'>#</a>
**`ClimateTools.frostdays`** &mdash; *Method*.



frostdays(data::Array, time::StepRange{Date,Base.Dates.Day})

FD, Number of frost days: Annual count of days when TN (daily minimum temperature) < 0 Celsius.

Let TN(i,j) be daily minimum temperature on day i in year j. Count the number of days where:

TN(i,j) < 0 Celsius.


<a target='_blank' href='https://github.com/Balinus/ClimateTools.jl/tree/64bde61fd1b4046c51b603edb114089b9e61186e/src/indices.jl#L42-L50' class='documenter-source'>source</a><br>

<a id='ClimateTools.icingdays-Tuple{Array{Float64,1},StepRange{Date,Base.Dates.Day}}' href='#ClimateTools.icingdays-Tuple{Array{Float64,1},StepRange{Date,Base.Dates.Day}}'>#</a>
**`ClimateTools.icingdays`** &mdash; *Method*.



icingdays(data::Array, time::StepRange{Date,Base.Dates.Day})

ID, Number of summer days: Annual count of days when TX (daily maximum temperature) < 0 degree Celsius.

Let TX(i,j) be daily maximum temperature on day i in year j. Count the number of days where:

TX(i,j) < 0 Celsius.


<a target='_blank' href='https://github.com/Balinus/ClimateTools.jl/tree/64bde61fd1b4046c51b603edb114089b9e61186e/src/indices.jl#L132-L140' class='documenter-source'>source</a><br>

<a id='ClimateTools.annualmin-Tuple{Array{Float64,1},StepRange{Date,Base.Dates.Day}}' href='#ClimateTools.annualmin-Tuple{Array{Float64,1},StepRange{Date,Base.Dates.Day}}'>#</a>
**`ClimateTools.annualmin`** &mdash; *Method*.



annualmin(data::Array, time::StepRange{Date,Base.Dates.Day})

AM, Value of annual minimum of array data.

Let data(i,j) be daily time serie on day i in year j. Extract the lowest value for year j.


<a target='_blank' href='https://github.com/Balinus/ClimateTools.jl/tree/64bde61fd1b4046c51b603edb114089b9e61186e/src/indices.jl#L332-L338' class='documenter-source'>source</a><br>

<a id='ClimateTools.annualmax-Tuple{Array{Float64,1},StepRange{Date,Base.Dates.Day}}' href='#ClimateTools.annualmax-Tuple{Array{Float64,1},StepRange{Date,Base.Dates.Day}}'>#</a>
**`ClimateTools.annualmax`** &mdash; *Method*.



annualmax(data::Array, time::StepRange{Date,Base.Dates.Day})

AM, Value of annual maximum of array data.

Let data(i,j) be daily time serie on day i in year j. Extract the highest value for year j.


<a target='_blank' href='https://github.com/Balinus/ClimateTools.jl/tree/64bde61fd1b4046c51b603edb114089b9e61186e/src/indices.jl#L289-L295' class='documenter-source'>source</a><br>

<a id='ClimateTools.tropicalnights-Tuple{Array{Float64,1},StepRange{Date,Base.Dates.Day}}' href='#ClimateTools.tropicalnights-Tuple{Array{Float64,1},StepRange{Date,Base.Dates.Day}}'>#</a>
**`ClimateTools.tropicalnights`** &mdash; *Method*.



tropicalnights(data::Array, time::StepRange{Date,Base.Dates.Day})

TropicalNights, Number of tropical nights: Annual count of days when TN (daily maximum temperature) > 20 degree Celsius.

Let TN(i,j) be daily minimum temperature on day i in year j. Count the number of days where:

TN(i,j) > 20 Celsius.


<a target='_blank' href='https://github.com/Balinus/ClimateTools.jl/tree/64bde61fd1b4046c51b603edb114089b9e61186e/src/indices.jl#L153-L161' class='documenter-source'>source</a><br>

<a id='ClimateTools.customthresover-Tuple{Array{Float64,1},StepRange{Date,Base.Dates.Day},Any}' href='#ClimateTools.customthresover-Tuple{Array{Float64,1},StepRange{Date,Base.Dates.Day},Any}'>#</a>
**`ClimateTools.customthresover`** &mdash; *Method*.



customthresover(data::Array, time::StepRange{Date,Base.Dates.Day}, thres)

customthresover, annual number of days over a specified threshold.

Let TS(i,j) be a daily time serie value on day i in year j. Count the number of days where:

TS(i,j) > thres.


<a target='_blank' href='https://github.com/Balinus/ClimateTools.jl/tree/64bde61fd1b4046c51b603edb114089b9e61186e/src/indices.jl#L198-L206' class='documenter-source'>source</a><br>

<a id='ClimateTools.customthresunder-Tuple{Array{Float64,1},StepRange{Date,Base.Dates.Day},Any}' href='#ClimateTools.customthresunder-Tuple{Array{Float64,1},StepRange{Date,Base.Dates.Day},Any}'>#</a>
**`ClimateTools.customthresunder`** &mdash; *Method*.



customthresunder(data::Array, time::StepRange{Date,Base.Dates.Day}, thres)

customthresover, annual number of days under a specified threshold.

Let TS(i,j) be a daily time serie value on day i in year j. Count the number of days where:

TS(i,j) < thres.


<a target='_blank' href='https://github.com/Balinus/ClimateTools.jl/tree/64bde61fd1b4046c51b603edb114089b9e61186e/src/indices.jl#L243-L251' class='documenter-source'>source</a><br>

<a id='ClimateTools.prcp1-Tuple{Array{Float64,1},StepRange{Date,Base.Dates.Day}}' href='#ClimateTools.prcp1-Tuple{Array{Float64,1},StepRange{Date,Base.Dates.Day}}'>#</a>
**`ClimateTools.prcp1`** &mdash; *Method*.



prcp1(data::Array, timevector::StepRange{Date,Base.Dates.Day})

Annual number with preciptation over 1 mm. This function returns a boolean vector. *true* if the data is higher or equal to 1 and *false* otherwise.


<a target='_blank' href='https://github.com/Balinus/ClimateTools.jl/tree/64bde61fd1b4046c51b603edb114089b9e61186e/src/indices.jl#L1-L5' class='documenter-source'>source</a><br>

<a id='ClimateTools.inpoly-Tuple{Any,Array{T,2}}' href='#ClimateTools.inpoly-Tuple{Any,Array{T,2}}'>#</a>
**`ClimateTools.inpoly`** &mdash; *Method*.



Determines if a point is inside a polygon.

  * p – point (x,y) or [x,y]
  * poly – polygon vertices [x1 x2 ... xn x1                           y1 y2 ... yn y1](a closed poly)

Returns true if point has an odd winding number.  This should label points as exterior which are inside outcrops.  See test for a test.


<a target='_blank' href='https://github.com/Balinus/ClimateTools.jl/tree/64bde61fd1b4046c51b603edb114089b9e61186e/src/functions.jl#L44-L54' class='documenter-source'>source</a><br>

<a id='ClimateTools.summerdays-Tuple{Array{Float64,1},StepRange{Date,Base.Dates.Day}}' href='#ClimateTools.summerdays-Tuple{Array{Float64,1},StepRange{Date,Base.Dates.Day}}'>#</a>
**`ClimateTools.summerdays`** &mdash; *Method*.



summerdays(data::Array, time::StepRange{Date,Base.Dates.Day})

SD, Number of summer days: Annual count of days when TX (daily maximum temperature) > 25 degree Celsius.

Let TX(i,j) be daily maximum temperature on day i in year j. Count the number of days where:

TX(i,j) >= 25 Celsius.


<a target='_blank' href='https://github.com/Balinus/ClimateTools.jl/tree/64bde61fd1b4046c51b603edb114089b9e61186e/src/indices.jl#L87-L95' class='documenter-source'>source</a><br>

<a id='ClimateTools.meshgrid-Tuple{Array{Float64,1},StepRange{Date,Base.Dates.Day}}' href='#ClimateTools.meshgrid-Tuple{Array{Float64,1},StepRange{Date,Base.Dates.Day}}'>#</a>
**`ClimateTools.meshgrid`** &mdash; *Method*.



This function creates a 2-D mesh-grid in a format consistent with Matlab's function meshgrid()

[X, Y] = meshgrid(XV, YV)

where XV and YV are vectors.


<a target='_blank' href='https://github.com/Balinus/ClimateTools.jl/tree/64bde61fd1b4046c51b603edb114089b9e61186e/src/functions.jl#L60-L66' class='documenter-source'>source</a><br>

<a id='ClimateTools.windnr-Tuple{Any,Array{T,2}}' href='#ClimateTools.windnr-Tuple{Any,Array{T,2}}'>#</a>
**`ClimateTools.windnr`** &mdash; *Method*.



Determines the winding number of a point and a polygon, i.e. how many times a polygon winds around the point.

It follows Dan Sunday: http://geomalgorithms.com/a03-_inclusion.html.


<a target='_blank' href='https://github.com/Balinus/ClimateTools.jl/tree/64bde61fd1b4046c51b603edb114089b9e61186e/src/functions.jl#L1-L6' class='documenter-source'>source</a><br>

<a id='ClimateTools.boxcar3-Tuple{AbstractArray}' href='#ClimateTools.boxcar3-Tuple{AbstractArray}'>#</a>
**`ClimateTools.boxcar3`** &mdash; *Method*.



This function creates a boxcar averager with a window length of 3

function boxcar3(A::AbstractArray)


<a target='_blank' href='https://github.com/Balinus/ClimateTools.jl/tree/64bde61fd1b4046c51b603edb114089b9e61186e/src/functions.jl#L86-L91' class='documenter-source'>source</a><br>


<a id='Index-1'></a>

## Index

- [`ClimateTools.annualmax`](index.md#ClimateTools.annualmax-Tuple{Array{Float64,1},StepRange{Date,Base.Dates.Day}})
- [`ClimateTools.annualmin`](index.md#ClimateTools.annualmin-Tuple{Array{Float64,1},StepRange{Date,Base.Dates.Day}})
- [`ClimateTools.boxcar3`](index.md#ClimateTools.boxcar3-Tuple{AbstractArray})
- [`ClimateTools.customthresover`](index.md#ClimateTools.customthresover-Tuple{Array{Float64,1},StepRange{Date,Base.Dates.Day},Any})
- [`ClimateTools.customthresunder`](index.md#ClimateTools.customthresunder-Tuple{Array{Float64,1},StepRange{Date,Base.Dates.Day},Any})
- [`ClimateTools.frostdays`](index.md#ClimateTools.frostdays-Tuple{Array{Float64,1},StepRange{Date,Base.Dates.Day}})
- [`ClimateTools.icingdays`](index.md#ClimateTools.icingdays-Tuple{Array{Float64,1},StepRange{Date,Base.Dates.Day}})
- [`ClimateTools.inpoly`](index.md#ClimateTools.inpoly-Tuple{Any,Array{T,2}})
- [`ClimateTools.meshgrid`](index.md#ClimateTools.meshgrid-Tuple{Array{Float64,1},StepRange{Date,Base.Dates.Day}})
- [`ClimateTools.prcp1`](index.md#ClimateTools.prcp1-Tuple{Array{Float64,1},StepRange{Date,Base.Dates.Day}})
- [`ClimateTools.summerdays`](index.md#ClimateTools.summerdays-Tuple{Array{Float64,1},StepRange{Date,Base.Dates.Day}})
- [`ClimateTools.tropicalnights`](index.md#ClimateTools.tropicalnights-Tuple{Array{Float64,1},StepRange{Date,Base.Dates.Day}})
- [`ClimateTools.windnr`](index.md#ClimateTools.windnr-Tuple{Any,Array{T,2}})

