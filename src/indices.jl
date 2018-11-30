# TODO add quantile estimation as indices
function buildarray_annual(C::ClimateTools.ClimGrid, dataout, numYears)
    lonsymbol = Symbol(C.dimension_dict["lon"])
    latsymbol = Symbol(C.dimension_dict["lat"])
    FD = AxisArray(dataout, Axis{lonsymbol}(C[1][Axis{lonsymbol}][:]), Axis{latsymbol}(C[1][Axis{latsymbol}][:]), Axis{:time}(Dates.Year.(numYears)))
    return FD
end

function buildarray_climato(C::ClimateTools.ClimGrid, dataout)
    lonsymbol = Symbol(C.dimension_dict["lon"])
    latsymbol = Symbol(C.dimension_dict["lat"])
    FD = AxisArray(dataout, Axis{lonsymbol}(C[1][Axis{lonsymbol}][:]), Axis{latsymbol}(C[1][Axis{latsymbol}][:]))
    return FD
end

function buildarray_orig(C::ClimateTools.ClimGrid, dataout)
    lonsymbol = Symbol(C.dimension_dict["lon"])
    latsymbol = Symbol(C.dimension_dict["lat"])
    FD = AxisArray(dataout, Axis{lonsymbol}(C[1][Axis{lonsymbol}][:]), Axis{latsymbol}(C[1][Axis{latsymbol}][:]), Axis{:time}(get_timevec(C)))
    return FD
end

function buildarray_resample(C::ClimateTools.ClimGrid, dataout, newtime)
    lonsymbol = Symbol(C.dimension_dict["lon"])
    latsymbol = Symbol(C.dimension_dict["lat"])
    FD = AxisArray(dataout, Axis{lonsymbol}(C[1][Axis{lonsymbol}][:]), Axis{latsymbol}(C[1][Axis{latsymbol}][:]), Axis{:time}(newtime))
    return FD
end

"""
    prcp1(C::ClimGrid)

Annual number with preciptation >= 1 mm. This function returns a ClimGrid.
"""
function prcp1(C::ClimGrid)
  @argcheck C[9] == "pr"
  years    = Dates.year.(C.data[Axis{:time}][:])
  numYears = unique(years)
  dataout  = zeros(Int64, (size(C.data, 1), size(C.data, 2), length(numYears)))
  datain   = C.data.data

  # Indice calculation
  Threads.@threads for i in 1:length(numYears)
    idx = searchsortedfirst(years, numYears[i]):searchsortedlast(years, numYears[i])
    Base.mapreducedim!(t -> t >= 1, +, view(dataout, :, :, i:i), view(datain, :,:, idx))
  end

  # Apply mask
  dataout = applymask(dataout, C.msk)

  # Build output AxisArray
  FD = buildarray_annual(C, dataout, numYears)

  # Return climGrid type containing the indice
  return ClimGrid(FD, longrid=C.longrid, latgrid=C.latgrid, msk=C.msk, grid_mapping=C.grid_mapping, dimension_dict=C.dimension_dict, timeattrib=C.timeattrib, model=C.model, frequency="year", experiment=C.experiment, run=C.run, project=C.project, institute=C.institute, filename=C.filename, dataunits="days", latunits=C.latunits, lonunits=C.lonunits, variable="prcp1", typeofvar=C.typeofvar, typeofcal=C.typeofcal, varattribs=C.varattribs, globalattribs=C.globalattribs)
end

"""
    frostdays(C::ClimGrid)

FD, Number of frost days: Annual count of days when TN (daily minimum temperature) < 0 Celsius.

Let TN[i,j] be daily minimum temperature on day i in year j. Count the number of days where:

  TN[i,j] < 0 Celsius.
"""
function frostdays(C::ClimGrid)
  @argcheck C[9] == "tasmin"
  years    = Dates.year.(C.data[Axis{:time}][:])
  numYears = unique(years)
  dataout  = zeros(Int64, (size(C.data, 1), size(C.data, 2), length(numYears)))
  datain   = C.data.data

  # Indice calculation
  Threads.@threads for i in 1:length(numYears)
    idx = searchsortedfirst(years, numYears[i]):searchsortedlast(years, numYears[i])
    Base.mapreducedim!(t -> t < 0, +, view(dataout, :, :, i:i), view(datain, :,:, idx))
  end

  # Apply mask
  dataout = applymask(dataout, C.msk)

  # Build output AxisArray
  FD = buildarray_annual(C, dataout, numYears)

  # Return climGrid type containing the indice
  return ClimGrid(FD, longrid=C.longrid, latgrid=C.latgrid, msk=C.msk, grid_mapping=C.grid_mapping, dimension_dict=C.dimension_dict, timeattrib=C.timeattrib, model=C.model, frequency="year", experiment=C.experiment, run=C.run, project=C.project, institute=C.institute, filename=C.filename, dataunits="days", latunits=C.latunits, lonunits=C.lonunits, variable="frostdays", typeofvar=C.typeofvar, typeofcal=C.typeofcal, varattribs=C.varattribs, globalattribs=C.globalattribs)
end

"""
    summerdays(C::ClimGrid)

SD, Number of summer days: Annual count of days when TX (daily maximum temperature) > 25 degree Celsius.

Let TX[i,j] be daily maximum temperature on day i in year j. Count the number of days where:

  TX[i,j] > 25 Celsius.
"""
function summerdays(C::ClimGrid)
  @argcheck C[9] == "tasmax"
  years    = Dates.year.(C.data[Axis{:time}][:])
  numYears = unique(years)
  dataout  = zeros(Int64, (size(C.data, 1), size(C.data, 2), length(numYears)))
  datain   = C.data.data

  # Indice calculation
  Threads.@threads for i in 1:length(numYears)
    idx = searchsortedfirst(years, numYears[i]):searchsortedlast(years, numYears[i])
    Base.mapreducedim!(t -> t > 25, +, view(dataout, :, :, i:i), view(datain, :,:, idx))
  end

  # Apply mask
  dataout = applymask(dataout, C.msk)

  # Build output AxisArray
  FD = buildarray_annual(C, dataout, numYears)

  # Return climGrid type containing the indice
  return ClimGrid(FD, longrid=C.longrid, latgrid=C.latgrid, msk=C.msk, grid_mapping=C.grid_mapping, dimension_dict=C.dimension_dict, timeattrib=C.timeattrib, model=C.model, frequency="year", experiment=C.experiment, run=C.run, project=C.project, institute=C.institute, filename=C.filename, dataunits="days", latunits=C.latunits, lonunits=C.lonunits, variable="summerdays", typeofvar=C.typeofvar, typeofcal=C.typeofcal, varattribs=C.varattribs, globalattribs=C.globalattribs)
end

"""
    icingdays(C::ClimGrid)

ID, Number of summer days: Annual count of days when TX (daily maximum temperature) < 0 degree Celsius.

Let TX[i,j] be daily maximum temperature on day i in year j. Count the number of days where:

  TX[i,j] < 0 Celsius.
"""
function icingdays(C::ClimGrid)
  @argcheck C[9] == "tasmax"
  years    = Dates.year.(C.data[Axis{:time}][:])
  numYears = unique(years)
  dataout  = zeros(Int64, (size(C.data, 1), size(C.data, 2), length(numYears)))
  datain   = C.data.data

  # Indice calculation
  Threads.@threads for i in 1:length(numYears)
    idx = searchsortedfirst(years, numYears[i]):searchsortedlast(years, numYears[i])
    Base.mapreducedim!(t -> t < 0, +, view(dataout, :, :, i:i), view(datain, :,:, idx))
  end

  # Apply mask
  dataout = applymask(dataout, C.msk)

  # Build output AxisArray
  FD = buildarray_annual(C, dataout, numYears)

  # Return climGrid type containing the indice
  return ClimGrid(FD, longrid=C.longrid, latgrid=C.latgrid, msk=C.msk, grid_mapping=C.grid_mapping, dimension_dict=C.dimension_dict, timeattrib=C.timeattrib, model=C.model, frequency="year", experiment=C.experiment, run=C.run, project=C.project, institute=C.institute, filename=C.filename, dataunits="days", latunits=C.latunits, lonunits=C.lonunits, variable="icingdays", typeofvar=C.typeofvar, typeofcal=C.typeofcal, varattribs=C.varattribs, globalattribs=C.globalattribs)
end


"""
    tropicalnights(C::ClimGrid)

TropicalNights, Number of tropical nights: Annual count of days when TN (daily maximum temperature) > 20 degree Celsius.

Let TN[i,j] be daily minimum temperature on day i in year j. Count the number of days where:

  TN[i,j] > 20 Celsius.
"""
function tropicalnights(C::ClimGrid)
  @argcheck C[9] == "tasmin"
  years    = Dates.year.(C.data[Axis{:time}][:])
  numYears = unique(years)
  dataout  = zeros(Int64, (size(C.data, 1), size(C.data, 2), length(numYears)))
  datain   = C.data.data

  # Indice calculation
  Threads.@threads for i in 1:length(numYears)
    idx = searchsortedfirst(years, numYears[i]):searchsortedlast(years, numYears[i])
    Base.mapreducedim!(t -> t > 20, +, view(dataout, :, :, i:i), view(datain, :,:, idx))
  end

  # Apply mask
  dataout = applymask(dataout, C.msk)

  # Build output AxisArray
  FD = buildarray_annual(C, dataout, numYears)

  # Return climGrid type containing the indice
  return ClimGrid(FD, longrid=C.longrid, latgrid=C.latgrid, msk=C.msk, grid_mapping=C.grid_mapping, dimension_dict=C.dimension_dict, timeattrib=C.timeattrib, model=C.model, frequency="year", experiment=C.experiment, run=C.run, project=C.project, institute=C.institute, filename=C.filename, dataunits="days", latunits=C.latunits, lonunits=C.lonunits, variable="tropicalnights", typeofvar=C.typeofvar, typeofcal=C.typeofcal, varattribs=C.varattribs, globalattribs=C.globalattribs)
end

"""
    customthresover(C::ClimGrid)

customthresover, annual number of days over a specified threshold.

Let TS[i,j] be a daily time serie value on day i in year j. Count the number of days where:

  TS[i,j] > thres.
"""
function customthresover(C::ClimGrid, thres)
  years    = Dates.year.(C.data[Axis{:time}][:])
  numYears = unique(years)
  dataout  = fill(NaN, (size(C.data, 1), size(C.data, 2), length(numYears)))
  datain   = C.data.data

  dataout_rshp = reshape(dataout, (size(dataout, 1)*size(dataout, 2), size(dataout, 3)))
  datain_rshp = reshape(datain, (size(datain, 1)*size(datain, 2), size(datain, 3)))

  # Indice calculation
  Threads.@threads for k in 1:size(datain_rshp, 1)
    for i in 1:length(numYears)
      idx = searchsortedfirst(years, numYears[i]):searchsortedlast(years, numYears[i])
      dataout_rshp[k, i] = sum(datain_rshp[k, idx] .> thres)
      # Base.mapreducedim!(t -> t > thres, +, view(dataout, :, :, i:i), view(datain, :,:, idx))
    end
  end

  # Apply mask
  dataout = applymask(dataout, C.msk)

  # Build output AxisArray
  FD = buildarray_annual(C, dataout, numYears)

  # Return climGrid type containing the indice
  return ClimGrid(FD, longrid=C.longrid, latgrid=C.latgrid, msk=C.msk, grid_mapping=C.grid_mapping, dimension_dict=C.dimension_dict, timeattrib=C.timeattrib, model=C.model, frequency="year", experiment=C.experiment, run=C.run, project=C.project, institute=C.institute, filename=C.filename, dataunits="days", latunits=C.latunits, lonunits=C.lonunits, variable=string("Days over ", thres, " ", C.dataunits), typeofvar=C.typeofvar, typeofcal=C.typeofcal, varattribs=C.varattribs, globalattribs=C.globalattribs)
end

"""
    customthresunder(C::ClimGrid)

customthresover, annual number of days under a specified threshold.

Let TS[i,j] be a daily time serie value on day i in year j. Count the number of days where:

  TS[i,j] < thres.
"""
function customthresunder(C::ClimGrid, thres)
  years    = Dates.year.(C.data[Axis{:time}][:])
  numYears = unique(years)
  dataout  = zeros(Int64, (size(C.data, 1), size(C.data, 2), length(numYears)))
  datain   = C.data.data

  # Indice calculation
  Threads.@threads for i in 1:length(numYears)
    idx = searchsortedfirst(years, numYears[i]):searchsortedlast(years, numYears[i])
    Base.mapreducedim!(t -> t < thres, +, view(dataout, :, :, i:i), view(datain, :,:, idx))
  end

  # Apply mask
  dataout = applymask(dataout, C.msk)

  # Build output AxisArray
  FD = buildarray_annual(C, dataout, numYears)

  # Return climGrid type containing the indice
  return ClimGrid(FD, longrid=C.longrid, latgrid=C.latgrid, msk=C.msk, grid_mapping=C.grid_mapping, dimension_dict=C.dimension_dict, timeattrib=C.timeattrib, model=C.model, frequency="year", experiment=C.experiment, run=C.run, project=C.project, institute=C.institute, filename=C.filename, dataunits="days", latunits=C.latunits, lonunits=C.lonunits, variable=string("Days under ", thres, " ", C.dataunits), typeofvar=C.typeofvar, typeofcal=C.typeofcal, varattribs=C.varattribs, globalattribs=C.globalattribs)
end

"""
    annualmax(C::ClimGrid)

Annual maximum of array data.

Let data[i,j] be daily time serie on day i in year j. Extract the highest value for year j.
"""
function annualmax(C::ClimGrid)
  years    = Dates.year.(C.data[Axis{:time}][:])
  numYears = unique(years)
  dataout  = zeros(Float64, (size(C.data, 1), size(C.data, 2), length(numYears)))
  datain   = C.data.data

  # Indice calculation
  Threads.@threads for i in 1:length(numYears)
    idx = searchsortedfirst(years, numYears[i]):searchsortedlast(years, numYears[i])
    Statistics.maximum!(view(dataout, :, :, i:i), view(datain, :,:, idx))
  end

  # Apply mask
  dataout = applymask(dataout, C.msk)

  # Build output AxisArray
  FD = buildarray_annual(C, dataout, numYears)

  # Return climGrid type containing the indice
  return ClimGrid(FD, longrid=C.longrid, latgrid=C.latgrid, msk=C.msk, grid_mapping=C.grid_mapping, dimension_dict=C.dimension_dict, timeattrib=C.timeattrib, model=C.model, frequency="year", experiment=C.experiment, run=C.run, project=C.project, institute=C.institute, filename=C.filename, dataunits=C.dataunits, latunits=C.latunits, lonunits=C.lonunits, variable="annualmax", typeofvar=C.typeofvar, typeofcal=C.typeofcal, varattribs=C.varattribs, globalattribs=C.globalattribs)
end


"""
    annualmin(C::ClimGrid)

Annual minimum of array data.

Let data[i,j] be daily time serie on day i in year j. Extract the lowest value for year j.
"""
function annualmin(C::ClimGrid)
  years    = Dates.year.(C.data[Axis{:time}][:])
  numYears = unique(years)
  dataout  = zeros(Float64, (size(C.data, 1), size(C.data, 2), length(numYears)))
  datain   = C.data.data

  # Indice calculation
  Threads.@threads for i in 1:length(numYears)
    idx = searchsortedfirst(years, numYears[i]):searchsortedlast(years, numYears[i])
    Base.minimum!(view(dataout, :, :, i:i), view(datain, :,:, idx))
  end

  # Apply mask
  dataout = applymask(dataout, C.msk)

  # Build output AxisArray
  FD = buildarray_annual(C, dataout, numYears)

  # Return climGrid type containing the indice
  return ClimGrid(FD, longrid=C.longrid, latgrid=C.latgrid, msk=C.msk, grid_mapping=C.grid_mapping, dimension_dict=C.dimension_dict, timeattrib=C.timeattrib, model=C.model, frequency="year", experiment=C.experiment, run=C.run, project=C.project, institute=C.institute, filename=C.filename, dataunits=C.dataunits, latunits=C.latunits, lonunits=C.lonunits, variable="annualmin", typeofvar=C.typeofvar, typeofcal=C.typeofcal, varattribs=C.varattribs, globalattribs=C.globalattribs)
end

"""
    annualmean(C::ClimGrid)

Annual mean of array data.

Let data[i,j] be daily time serie on day i in year j. Calculate the mean value for year j.
"""
function annualmean(C::ClimGrid)
  years    = Dates.year.(C.data[Axis{:time}][:])
  numYears = unique(years)
  dataout  = zeros(Float64, (size(C.data, 1), size(C.data, 2), length(numYears)))
  datain   = C.data.data

  # Indice calculation
  Threads.@threads for i in 1:length(numYears)
    idx = searchsortedfirst(years, numYears[i]):searchsortedlast(years, numYears[i])
    Statistics.mean!(view(dataout, :, :, i:i), view(datain, :,:, idx))
  end

  # Apply mask
  dataout = applymask(dataout, C.msk)

  # Build output AxisArray
  FD = buildarray_annual(C, dataout, numYears)

  # Return climGrid type containing the indice
  return ClimGrid(FD, longrid=C.longrid, latgrid=C.latgrid, msk=C.msk, grid_mapping=C.grid_mapping, dimension_dict=C.dimension_dict, timeattrib=C.timeattrib, model=C.model, frequency="year", experiment=C.experiment, run=C.run, project=C.project, institute=C.institute, filename=C.filename, dataunits=C.dataunits, latunits=C.latunits, lonunits=C.lonunits, variable="annualmean", typeofvar=C.typeofvar, typeofcal=C.typeofcal, varattribs=C.varattribs, globalattribs=C.globalattribs)
end

"""
    annualsum(C::ClimGrid)

Annual sum of array data.

Let data[i,j] be daily time serie on day i in year j. Sums daily values for year j.
"""
function annualsum(C::ClimGrid)
  years    = Dates.year.(C.data[Axis{:time}][:])
  numYears = unique(years)
  dataout  = zeros(Float64, (size(C.data, 1), size(C.data, 2), length(numYears)))
  datain   = C.data.data

  # Indice calculation
  Threads.@threads for i in 1:length(numYears)
    idx = searchsortedfirst(years, numYears[i]):searchsortedlast(years, numYears[i])
    Base.sum!(view(dataout, :, :, i:i), view(datain, :,:, idx))
  end

  # Apply mask
  dataout = applymask(dataout, C.msk)

  # Build output AxisArray
  FD = buildarray_annual(C, dataout, numYears)

  # Return climGrid type containing the indice
  return ClimGrid(FD, longrid=C.longrid, latgrid=C.latgrid, msk=C.msk, grid_mapping=C.grid_mapping, dimension_dict=C.dimension_dict, timeattrib=C.timeattrib, model=C.model, frequency="year", experiment=C.experiment, run=C.run, project=C.project, institute=C.institute, filename=C.filename, dataunits=C.dataunits, latunits=C.latunits, lonunits=C.lonunits, variable="annualsum", typeofvar=C.typeofvar, typeofcal=C.typeofcal, varattribs=C.varattribs, globalattribs=C.globalattribs)
end


"""
    spei()

Returns the The Standardized Precipitation Evapotranspiration Index (SPEI).

**Reference**
Vicente-Serrano, S. M., Beguería, S., & López-Moreno, J. I. (2010). A Multiscalar Drought Index Sensitive to Global Warming: The Standardized Precipitation Evapotranspiration Index. Journal of Climate, 23(7), 1696–1718. https://doi.org/10.1175/2009JCLI2909.1
"""

function spei()

end
