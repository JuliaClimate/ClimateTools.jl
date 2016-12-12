
<a id='Functions-in-ClimateTools-package-1'></a>

# Functions in ClimateTools package


Description of the functions


<a id='Functions-1'></a>

## Functions

<a id='ClimateTools.prcp1' href='#ClimateTools.prcp1'>#</a>
**`ClimateTools.prcp1`** &mdash; *Function*.



prcp1(data::Array, timevector::StepRange{Date,Base.Dates.Day})

Annual number with preciptation over 1 mm. This function returns a boolean vector. *true* if the data is higher or equal to 1 and *false* otherwise.


<a target='_blank' href='https://github.com/Balinus/ClimateTools.jl/tree/9fcfec0ee35bf3383d91d36a3583abbecf0c52b1/src/indices.jl#L1-L5' class='documenter-source'>source</a><br>

<a id='ClimateTools.inpoly' href='#ClimateTools.inpoly'>#</a>
**`ClimateTools.inpoly`** &mdash; *Function*.



Determines if a point is inside a polygon.

  * p – point (x,y) or [x,y]
  * poly – polygon vertices [x1 x2 ... xn x1                           y1 y2 ... yn y1](a closed poly)

Returns true if point has an odd winding number.  This should label points as exterior which are inside outcrops.  See test for a test.


<a target='_blank' href='https://github.com/Balinus/ClimateTools.jl/tree/9fcfec0ee35bf3383d91d36a3583abbecf0c52b1/src/functions.jl#L44-L54' class='documenter-source'>source</a><br>

<a id='ClimateTools.windnr' href='#ClimateTools.windnr'>#</a>
**`ClimateTools.windnr`** &mdash; *Function*.



Determines the winding number of a point and a polygon, i.e. how many times a polygon winds around the point.

It follows Dan Sunday: http://geomalgorithms.com/a03-_inclusion.html.


<a target='_blank' href='https://github.com/Balinus/ClimateTools.jl/tree/9fcfec0ee35bf3383d91d36a3583abbecf0c52b1/src/functions.jl#L1-L6' class='documenter-source'>source</a><br>

<a id='ClimateTools.boxcar3' href='#ClimateTools.boxcar3'>#</a>
**`ClimateTools.boxcar3`** &mdash; *Function*.



This function creates a boxcar averager with a window length of 3

function boxcar3(A::AbstractArray)


<a target='_blank' href='https://github.com/Balinus/ClimateTools.jl/tree/9fcfec0ee35bf3383d91d36a3583abbecf0c52b1/src/functions.jl#L83-L88' class='documenter-source'>source</a><br>

<a id='ClimateTools.meshgrid' href='#ClimateTools.meshgrid'>#</a>
**`ClimateTools.meshgrid`** &mdash; *Function*.



This function creates a 2-D mesh-grid in a format consistent with Matlab's function meshgrid()

[X, Y] = meshgrid(XV, YV)

where XV and YV are vectors.


<a target='_blank' href='https://github.com/Balinus/ClimateTools.jl/tree/9fcfec0ee35bf3383d91d36a3583abbecf0c52b1/src/functions.jl#L57-L63' class='documenter-source'>source</a><br>

