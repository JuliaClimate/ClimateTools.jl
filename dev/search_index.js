var documenterSearchIndex = {"docs": [

{
    "location": "#",
    "page": "Home",
    "title": "Home",
    "category": "page",
    "text": ""
},

{
    "location": "#Home-1",
    "page": "Home",
    "title": "Home",
    "category": "section",
    "text": ""
},

{
    "location": "#Overview-1",
    "page": "Home",
    "title": "Overview",
    "category": "section",
    "text": "ClimateTools.jl is a collection of commonly-used tools in Climate science. Basics of climate field analysis will be covered, with some forays into exploratory techniques. The package is aimed to ease the typical steps of analysis climate models outputs from netCDF files that follows Climate Forecast conventions and the creation of climate scenarios.The package is registered on METADATA.jl and can be added with] add ClimateTools\nusing ClimateTools"
},

{
    "location": "#Philosophy-1",
    "page": "Home",
    "title": "Philosophy",
    "category": "section",
    "text": "The idea behind ClimateTools is that most, if not all, climate fields can be represented by a 2D (e.g. topography), 3D (e.g. air temperature) or 4D (e.g. winds at multiple levels) grids that are georeferenced. Those grids are named ClimGrid in ClimateTools. Every functions acts on such structure and returns a similar structure. The ClimGrid structure contains all elements needed to be manipulated: latitude, longitude, calendars, variable attributes, etc. that was either available in the original netCDF file or that was inferred by the metadata. Note that a ClimGrid is defined for a single variable.The metadata follows the various transformations and is modified when necessary. For example, calculating the annual number of days with precipitation higher than 1mm will modify the variable name from pr (for precipitation) to prcp1, the name of the indicator. It will not, however, modify the base variable type (it will remain pr)."
},

{
    "location": "#Notes-1",
    "page": "Home",
    "title": "Notes",
    "category": "section",
    "text": "Where possible, functions are coded to use multiple threads. To gain maximum performance, use (bash shell) export JULIA_NUM_THREADS=n, where n is the number of threads. To get an idea of the number of threads you can use type (in Julia) Sys.THREADS. This is especially useful for climate indices, bias correction and regridding."
},

{
    "location": "#Features-1",
    "page": "Home",
    "title": "Features",
    "category": "section",
    "text": "Climate scenarios creation\nExtraction and visualization of CF-compliant netCDF datasets\nCustom user-provided polygons and start and end date for localized studies\nClimate indices from The joint CCl/CLIVAR/JCOMM Expert Team (ET) on Climate Change Detection and Indices (ETCCDI) as well as custom climate indices\nRegridding of a datasets onto another grid\nPost-processing of climate timeseries using Quantile-Quantile mapping method (cf. Themeßl et al. 2012, Piani et al. 2010)\nExportation of results to a CF-compliant netCDF file\nSupport for typical climate models calendars: 360day, 365day, Standard, Prolectip Gregorian through NCDatasets.jl.\nSupport for physical units through the Unitful.jl package."
},

{
    "location": "#Contributors-1",
    "page": "Home",
    "title": "Contributors",
    "category": "section",
    "text": "If you\'d like to have other climate indices coded, please, submit them through a Pull Request! I\'d be more than happy to include them. Alternatively, provide the equation in Issues."
},

{
    "location": "#TO-DO-1",
    "page": "Home",
    "title": "TO-DO",
    "category": "section",
    "text": "Dashboard tool. This will return the main characteristics of a ClimGrid: maps of minimum, maximum and mean climatological values, seasonal cycle, timeseries of annual maximum, minimum and mean values, etc...\nExtreme value theory analysis"
},

{
    "location": "#Documentation-1",
    "page": "Home",
    "title": "Documentation",
    "category": "section",
    "text": "This documentation was built using Documenter.jl.using Dates # hide\nprintln(\"Documentation built $(Dates.now()) with Julia $(VERSION)\") # hide"
},

{
    "location": "installation/#",
    "page": "Installation",
    "title": "Installation",
    "category": "page",
    "text": ""
},

{
    "location": "installation/#Installation-1",
    "page": "Installation",
    "title": "Installation",
    "category": "section",
    "text": ""
},

{
    "location": "installation/#Approach-no.-1-Use-main-system-python-distribution-1",
    "page": "Installation",
    "title": "Approach no. 1 Use main system python distribution",
    "category": "section",
    "text": "ClimateTools need some Python dependencies for mapping purpose. To ensure that ClimateTools works properly, it is recommended to use a Python distribution that can properly load the following python modules and build PyCall with the same python distribution.1.1 Python dependenciesmatplotlib (tested with version 2.0.1)\nbasemap (tested with version 1.0.7)\nscipy (tested with version 1.0.1)\ncmoceanNote. Installing Basemap for python 3.6+ seems problematic.1.2 Building PyCall After the confirmation that the Python dependencies can be loaded in Python, the user needs to build PyCall with the same Python version. Alternatively, if PyCall is already built, it may be only a matter of installing the Python dependencies with the PyCall\'s Python version by using pip.ENV[\"PYTHON\"]=\"path_to_python_distribution\"\npkg> build PyCall"
},

{
    "location": "installation/#Approach-no.-2.-Build-a-python-virtual-environment-and-link-PyCall.jl-to-it-1",
    "page": "Installation",
    "title": "Approach no. 2. Build a python virtual environment and link PyCall.jl to it",
    "category": "section",
    "text": "One approach to ensure that the right python dependencies are installed is to use a virtual environment. More information can be found in PyCall documentation.2.1 Create a virtual environment with Python 2.7.x.$ virtualenv --python=/usr/bin/python2 /path/to/venv\n$ /path/to/venv/bin/python -m pip install numpy\n$ /path/to/venv/bin/python -m pip install scipy\n$ /path/to/venv/bin/python -m pip install matplotlib\n$ /path/to/venv/bin/python -m pip install https://github.com/matplotlib/basemap/archive/v1.0.7rel.tar.gz\n$ /path/to/venv/bin/python -m pip install git+https://github.com/matplotlib/cmocean2.2 Testing Python installationLaunch newly created virtual env python#bash\n$ /path/to/venv/bin/python # launch virtual env pythonEnsure that you can load the appropriate python packages inside the python interpreter.#python\n>>> import mpl_toolkits.basemap as basemap\n>>> import matplotlib.pyplot as plt\n>>> import cmocean as cm\n>>> import scipy as sc2.3 Build PyCall with the new venv pythonOnce the virtual environment python is tested, it\'s a matter of telling PyCall to use this distribution.# julia\njulia> ENV[\"PYTHON\"] = \"/path/to/venv/bin/python\"\njulia> using Pkg;Pkg.build(\"PyCall\")\njulia> exit()\n# re-enter julia\njulia> using ClimateTools\njulia> using Pkg; Pkg.test(\"ClimateTools\")"
},

{
    "location": "installation/#Installing-ClimateTools.jl-1",
    "page": "Installation",
    "title": "Installing ClimateTools.jl",
    "category": "section",
    "text": "pkg> add ClimateTools # Tagged release"
},

{
    "location": "loadingfile/#",
    "page": "-",
    "title": "-",
    "category": "page",
    "text": ""
},

{
    "location": "loadingfile/#Reading-a-NetCDF-file-1",
    "page": "-",
    "title": "Reading a NetCDF file",
    "category": "section",
    "text": "The entry point of ClimateTools is to load data with the load function. The return structure of the load function is a in-memory representation of the variable contained in the netCDF file.C = load(filename::String, vari::String; poly::Array, data_units::String, start_date::Tuple, end_date::Tuple, dimension::Bool=true)load return a ClimGrid type. The ClimGrid represent a single variable. By default, the function tries to attach physical units to the data array by using the Unitful.jl package. The advantage behind physical units is that one can subtract a ClimGrid with Kelvin unit with a ClimGrid with Celsius unit and get coherent results. Be warned that some operations on some units are not allowed (you cannot \"add\" Celsius for instance). In the event that a user wants to do some calculations without physical logic, it is possible to load the dataset without the units by specifying dimension=false argument.Using the optional poly argument, the user can provide a polygon and the returned ClimGrid will only contains the grid points inside the provided polygon. The polygon provided should be in the -180, +180 longitude format. If the polygon crosses the International Date Line, the polygon should be splitted in multiple parts (i.e. multi-polygons).start_date and end_date can also be provided. It is useful when climate simulations file spans multiple decades/centuries and only a temporal subset is needed. Dates should be provided as a Tuple of the form (year, month, day, hour, minute, seconds), where only year is mandatory (e.g. (2000,) can be provided and will defaults to (2000, 01, 01)).For some variable, the optional keyword argument data_units can be provided. For example, precipitation in climate models are usually provided as kg/m^2/s. By specifying data_units = mm, the load function returns accumulation at the data time resolution. Similarly, the user can provide Celsius as data_units and load will return Celsius instead of Kelvin.struct ClimGrid\n  data::AxisArray # Data\n  longrid::AbstractArray{N,2} where N # the longitude grid\n  latgrid::AbstractArray{N,2} where N # the latitude grid\n  msk::Array{N, 2} where N # Data mask (NaNs and 1.0)\n  grid_mapping::Dict#{String, Any} # bindings for native grid\n  dimension_dict::Dict\n  model::String\n  frequency::String # Day, month, years\n  experiment::String # Historical, RCP4.5, RCP8.5, etc.\n  run::String\n  project::String # CORDEX, CMIP5, etc.\n  institute::String # UQAM, DMI, etc.\n  filename::String # Path of the original file\n  dataunits::String # Celsius, kelvin, etc.\n  latunits::String # latitude coordinate unit\n  lonunits::String # longitude coordinate unit\n  variable::String # Type of variable (i.e. can be the same as \"typeofvar\", but it is changed when calculating indices)\n  typeofvar::String # Variable type (e.g. tasmax, tasmin, pr)\n  typeofcal::String # Calendar type\n  timeattrib::Dict # Time attributes (e.g. days since ... )\n  varattribs::Dict # Variable attributes dictionary\n  globalattribs::Dict # Global attributes dictionary\nend"
},

{
    "location": "loadingfile/#Subsetting-1",
    "page": "-",
    "title": "Subsetting",
    "category": "section",
    "text": "Once the data is loaded in a ClimGrid struct, options to further subset the data are available."
},

{
    "location": "loadingfile/#Spatial-1",
    "page": "-",
    "title": "Spatial",
    "category": "section",
    "text": "spatialsubset function acts on ClimGrid type and subset the data through a spatial subset using a provided polygon. The function returns a ClimGrid. Polygons needs to be on a -180, +180 longitude coordinates, as data coordinates defaults to such grid. For instance, global models are often on a 0-360 degrees grid but the load function shift the data onto a -180,+180 coordinates.C = spatialsubset(C::ClimGrid, poly:Array{N, 2} where N)"
},

{
    "location": "loadingfile/#Temporal-1",
    "page": "-",
    "title": "Temporal",
    "category": "section",
    "text": "Temporal subset of the data is also possible with the temporalsubset function:C = temporalsubset(C::ClimGrid, startdate::Tuple, enddate::Tuple)"
},

{
    "location": "loadingfile/#Discontinuous-temporal-(e.g.-resampling)-1",
    "page": "-",
    "title": "Discontinuous temporal (e.g. resampling)",
    "category": "section",
    "text": "It is also possible to only keep a given non-continuous period for a given timeframe. For example, we might be interested in keeping only northern summer months (June-July-August) from a continuous ClimGrid covering 1961-2100. resample returns such a subsetted ClimGrid.Csub = resample(C, \"JJA\") # hardcoded ClimateTools\'s season\nCsub = resample(C, 6, 8) # custom subset example for June-July-August\nCsub = resample(C, 1, 2) # custom subset example for January-February"
},

{
    "location": "indices/#",
    "page": "Climate Indices",
    "title": "Climate Indices",
    "category": "page",
    "text": ""
},

{
    "location": "indices/#Climate-Indices-1",
    "page": "Climate Indices",
    "title": "Climate Indices",
    "category": "section",
    "text": ""
},

{
    "location": "indices/#Indices-1",
    "page": "Climate Indices",
    "title": "Indices",
    "category": "section",
    "text": "More than 20 climate indices are available in the package, such as the annual number of tropical nights, annual maximum and minimum, etc. You can calculate such indices simply with:ind = annualmax(C::ClimGrid)Which returns another ClimGrid. You can also map this ClimGrid with the mapclimgrid function and returns the climatological mean of the annual maximum (e.g. daily precipitation in the example below). A list of indices can be found in the documentation and in the functions.jl source code.mapclimgrid(ind) # mapping the indice previously calculated(Image: BNU-ESM)"
},

{
    "location": "indices/#Ensemble-mean-1",
    "page": "Climate Indices",
    "title": "Ensemble mean",
    "category": "section",
    "text": "You can calculate the ensemble mean with ensemble_mean function, where the input argument is an array of ClimGrids.Abstract example:C_model1 = ClimGrid(...) # model #1\nC_model2 = ClimGrid(...) # model #2\nens = [C_model1, C_model2] # Create an Array of ClimGrids\nE = ensemble_mean(ens) # Returns the mean of all models climatologies"
},

{
    "location": "indices/#ClimateTools.annualmax",
    "page": "Climate Indices",
    "title": "ClimateTools.annualmax",
    "category": "function",
    "text": "annualmax(C::ClimGrid)\n\nAnnual maximum of array data.\n\nLet data[i,j] be daily time serie on day i in year j. Extract the highest value for year j.\n\n\n\n\n\n"
},

{
    "location": "indices/#ClimateTools.annualmean",
    "page": "Climate Indices",
    "title": "ClimateTools.annualmean",
    "category": "function",
    "text": "annualmean(C::ClimGrid)\n\nAnnual mean of array data.\n\nLet data[i,j] be daily time serie on day i in year j. Calculate the mean value for year j.\n\n\n\n\n\n"
},

{
    "location": "indices/#ClimateTools.annualmin",
    "page": "Climate Indices",
    "title": "ClimateTools.annualmin",
    "category": "function",
    "text": "annualmin(C::ClimGrid)\n\nAnnual minimum of array data.\n\nLet data[i,j] be daily time serie on day i in year j. Extract the lowest value for year j.\n\n\n\n\n\n"
},

{
    "location": "indices/#ClimateTools.annualsum",
    "page": "Climate Indices",
    "title": "ClimateTools.annualsum",
    "category": "function",
    "text": "annualsum(C::ClimGrid)\n\nAnnual sum of array data.\n\nLet data[i,j] be daily time serie on day i in year j. Sums daily values for year j.\n\n\n\n\n\n"
},

{
    "location": "indices/#ClimateTools.approx_surfacepressure",
    "page": "Climate Indices",
    "title": "ClimateTools.approx_surfacepressure",
    "category": "function",
    "text": "  approx_surfacepressure(sealevel_pressure::ClimGrid, orography::ClimGrid, daily_temperature::ClimGrid)\n\nReturns the approximated surface pressure (sp) (Pa) using sea level pressure (psl) (Pa), orography (orog) (m), and daily mean temperature (tas) (K).\n\nsp = psl * 10^x\n\nwhere x = frac-orog18400 * tas  27315\n\n\n\n\n\n"
},

{
    "location": "indices/#ClimateTools.customthresover",
    "page": "Climate Indices",
    "title": "ClimateTools.customthresover",
    "category": "function",
    "text": "customthresover(C::ClimGrid, thres)\n\ncustomthresover, annual number of days over a specified threshold.\n\nLet TS[i,j] be a daily time serie value on day i in year j. Count the number of days where:\n\nTS[i,j] > thres.\n\nNote. The threshold needs to have units specified. For example:\n\njulia> using Unitful: @u_str, °C julia> thres = 15u\"°C\" 15 °C\n\n\n\n\n\n"
},

{
    "location": "indices/#ClimateTools.customthresunder",
    "page": "Climate Indices",
    "title": "ClimateTools.customthresunder",
    "category": "function",
    "text": "customthresunder(C::ClimGrid, thres::Quantity)\n\ncustomthresunder, annual number of days under a specified threshold.\n\nLet TS[i,j] be a daily time serie value on day i in year j. Count the number of days where:\n\nTS[i,j] < thres.\n\nNote. The threshold needs to have units specified. For example:\n\njulia> using Unitful: @u_str, °C julia> thres = 15u\"°C\" 15 °C\n\n\n\n\n\n"
},

{
    "location": "indices/#ClimateTools.diurnaltemperature",
    "page": "Climate Indices",
    "title": "ClimateTools.diurnaltemperature",
    "category": "function",
    "text": "diurnaltemperature(temperatureminimum::ClimGrid, temperaturemaximum::ClimGrid, α::Float64)\n\nReturns an estimation of the diurnal temperature (temperature between 7:00 (7am) and 17:00 (5pm)). The estimation is a linear combination of the daily minimum temperature (temperatureminimum) and daily maximum temperature (temperaturemaximum). The value of α has to be estimated seperatly from observations and depends on the location. The daily max and min must be in the same unit and in Celsius or Kelvin The diurnal temperature returned is in the same units as the daily minimum temperature and daily maximum temperature.\n\nTdiu = α * Tmin + (1 - α) * Tmax\n\n\n\n\n\n"
},

{
    "location": "indices/#ClimateTools.icingdays",
    "page": "Climate Indices",
    "title": "ClimateTools.icingdays",
    "category": "function",
    "text": "icingdays(C::ClimGrid)\n\nID, Number of summer days: Annual count of days when TX (daily maximum temperature) < 0 degree Celsius.\n\nLet TX[i,j] be daily maximum temperature on day i in year j. Count the number of days where:\n\nTX[i,j] < 0 Celsius.\n\n\n\n\n\n"
},

{
    "location": "indices/#ClimateTools.frostdays",
    "page": "Climate Indices",
    "title": "ClimateTools.frostdays",
    "category": "function",
    "text": "frostdays(C::ClimGrid)\n\nFD, Number of frost days: Annual count of days when TN (daily minimum temperature) < 0 Celsius.\n\nLet TN[i,j] be daily minimum temperature on day i in year j. Count the number of days where:\n\nTN[i,j] < 0 Celsius.\n\n\n\n\n\n"
},

{
    "location": "indices/#ClimateTools.prcp1",
    "page": "Climate Indices",
    "title": "ClimateTools.prcp1",
    "category": "function",
    "text": "prcp1(C::ClimGrid)\n\nAnnual number with preciptation >= 1 mm. This function returns a ClimGrid. Input data should be in mm.\n\n\n\n\n\n"
},

{
    "location": "indices/#ClimateTools.summerdays",
    "page": "Climate Indices",
    "title": "ClimateTools.summerdays",
    "category": "function",
    "text": "summerdays(C::ClimGrid)\n\nSD, Number of summer days: Annual count of days when TX (daily maximum temperature) > 25 degree Celsius.\n\nLet TX[i,j] be daily maximum temperature on day i in year j. Count the number of days where:\n\nTX[i,j] > 25 Celsius.\n\n\n\n\n\n"
},

{
    "location": "indices/#ClimateTools.tropicalnights",
    "page": "Climate Indices",
    "title": "ClimateTools.tropicalnights",
    "category": "function",
    "text": "tropicalnights(C::ClimGrid)\n\nTropicalNights, Number of tropical nights: Annual count of days when TN (daily maximum temperature) > 20 degree Celsius.\n\nLet TN[i,j] be daily minimum temperature on day i in year j. Count the number of days where:\n\nTN[i,j] > 20 Celsius.\n\n\n\n\n\n"
},

{
    "location": "indices/#ClimateTools.vaporpressure",
    "page": "Climate Indices",
    "title": "ClimateTools.vaporpressure",
    "category": "function",
    "text": "vaporpressure(surface_pressure::ClimGrid, specific_humidity::ClimGrid)\n\nReturns the vapor pressure (vp) (Pa) based on the surface pressure (sp) (Pa) and the specific humidity (q).\n\nvp = fracq * spq+0622\n\n\n\n\n\nvaporpressure(specific_humidity::ClimGrid, sealevel_pressure::ClimGrid, orography::ClimGrid, daily_temperature::ClimGrid)\n\nReturns the vapor pressure (vp) (Pa) estimated with the specific humidity (q), the sea level pressure (psl) (Pa), the orography (orog) (m) and the daily mean temperature (tas) (K). An approximation of the surface pressure is first computed by using the sea level pressure, orography and the daily mean temperature (see approx_surfacepressure). Then, vapor pressure is calculated by:\n\nvp = fracq * spq+0622\n\n\n\n\n\n"
},

{
    "location": "indices/#ClimateTools.wbgt",
    "page": "Climate Indices",
    "title": "ClimateTools.wbgt",
    "category": "function",
    "text": "wbgt(diurnal_temperature::ClimGrid, vapor_pressure::ClimGrid)\n\nReturns the simplified wet-bulb global temperature (wbgt) (Celsius) calculated using the vapor pressure (Pa) of the day and the estimated mean diurnal temperature (Celsius; temperature between 7:00 (7am) and 17:00 (5pm)).\n\nwbgt = 0567 * Tday + 000393 * vp + 394\n\n\n\n\n\n"
},

{
    "location": "indices/#Climate-Indices-2",
    "page": "Climate Indices",
    "title": "Climate Indices",
    "category": "section",
    "text": "Here\'s a list of climate indices currently provided by ClimateTools. This list may not be always up-to-date. See here for all exported functions.annualmax\nannualmean\nannualmin\nannualsum\napprox_surfacepressure\ncustomthresover\ncustomthresunder\ndiurnaltemperature\nicingdays\nfrostdays\nprcp1\nsummerdays\ntropicalnights\nvaporpressure\nwbgt"
},

{
    "location": "interpolation/#",
    "page": "Interpolation",
    "title": "Interpolation",
    "category": "page",
    "text": ""
},

{
    "location": "interpolation/#Interpolation-1",
    "page": "Interpolation",
    "title": "Interpolation",
    "category": "section",
    "text": "A typical step in climate analysis is to interpolate a given grid onto another grid. ClimateTools provides such a tool by wrapping Scipy griddata function. It is intended for visualization or as a 1st step before bias-correcting the ClimGrid dataset.regrid function will interpolate the data contained in ClimGrid A into the coordinates of ClimGrid B and returns a new ClimGrid C which contains the interpolated data of A into the grid of B.C = regrid(A::ClimGrid, B::ClimGrid)It is also possible to interpolate a ClimGrid onto specified longitude and latitude vectors and arrays.C = regrid(A::ClimGrid, lon::AbstractArray{N, T} where N where T, lat::AbstractArray{N, T} where N where T; dimx=[], dimy=[], method::String=\"linear\", min=[], max=[])In the case a longitude and latitude 2D array is provided, the user needs to provide the dimension vectors for x and y."
},

{
    "location": "biascorrection/#",
    "page": "Bias correction",
    "title": "Bias correction",
    "category": "page",
    "text": ""
},

{
    "location": "biascorrection/#Bias-correction-1",
    "page": "Bias correction",
    "title": "Bias correction",
    "category": "section",
    "text": "Quantile-quantile mapping (Themeßl et al. 2012, Grenier et al. 2015) is provided with ClimateTools.jl through the function qqmap.qqmap(obs::ClimGrid, ref::ClimGrid, fut::ClimGrid; method::String=\"Additive\", detrend::Bool=true, window::Int=15, rankn::Int=50, thresnan::Float64=0.1, keep_original::Bool=false, interp = Linear(), extrap = Flat())More information can be found in these references.Themeßl, Matthias Jakob, Andreas Gobiet, and Georg Heinrich. 2012. “Empirical-Statistical Downscaling and Error Correction of Regional Climate Models and Its Impact on the Climate Change Signal.” Climatic Change 112 (2). Springer: 449–68.Grenier, Patrick, Ramón de Elía, and Diane Chaumont. 2015. “Chances of Short-Term Cooling Estimated from a Selection of CMIP5-Based Climate Scenarios during 2006-2035 over Canada.” Journal of Climate, January 2015. American Meteorological Society. doi:10.1175/JCLI-D-14-00224.1."
},

{
    "location": "maps/#",
    "page": "Visualization",
    "title": "Visualization",
    "category": "page",
    "text": "CurrentModule = ClimateTools"
},

{
    "location": "maps/#Visualization-1",
    "page": "Visualization",
    "title": "Visualization",
    "category": "section",
    "text": ""
},

{
    "location": "maps/#ClimateTools.mapclimgrid",
    "page": "Visualization",
    "title": "ClimateTools.mapclimgrid",
    "category": "function",
    "text": "mapclimgrid(C::ClimGrid; region::String=\"auto\", poly, level, mask, caxis, start_date::Tuple, end_date::Tuple, titlestr::String, surface::Symbol, cm::String=\"\", ncolors::Int, center_cs::Bool, filename::String, cs_label::String)\n\nMaps the time-mean average of ClimGrid C. If a filename is provided, the figure is saved in a png format.\n\nOptional keyworkd includes precribed regions (keyword region, see list below), spatial clipping by polygon (keyword poly) or mask (keyword mask, an array of NaNs and 1.0 of the same dimension as the data in ClimGrid C), startdate and enddate. For 4D data, keyword level is used to map a given level (defaults to 1). caxis is used to limit the colorscale. cm is used to manually set the colorscale (see Python documentation for native colorscale keyword), ncolors is used to set the number of color classes (defaults to 12). Set center_cs to true to center the colorscale (useful for divergent results, such as anomalies, positive/negative temprature). cs_label is used for custom colorscale label.\n\nArguments for keyword region (and shortcuts)\n\nEurope (\"EU\")\nNorthAmerica (\"NA\")\nCanada (\"CA\")\nQuebec, QuebecNSP (\"QC\", \"QCNSP\")\nAmericas (\"Ams\")\nWorld, WorldAz, WorldEck4 (\"W\", \"Waz\", \"Weck4\")\nGreenwich (\"Gr\")\n\nArguments for keyword surface\n\n:contour\n:contourf\n:pcolormesh\n\n\n\n\n\nmapclimgrid(; region::String=\"auto\", poly, level, mask, caxis, start_date::Date, end_date::Date)\n\nEmpty map generator, when called without a ClimGrid as the positional argument.\n\n\n\n\n\n"
},

{
    "location": "maps/#Maps-1",
    "page": "Visualization",
    "title": "Maps",
    "category": "section",
    "text": "Mapping a ClimGrid is done by using the mapclimgrid function.C = load(filenc, \"pr\", data_units=\"mm\")\nmapclimgrid(C)(Image: CanESM2)mapclimgridNote that the function plots the climatological mean of the provided ClimGrid. Multiple options are available for region: World, Canada, Quebec, WorldAz, WorldEck4, ..., and the default auto which use the maximum and minimum of the lat-long coordinates inside the ClimGrid structure (see the documentation of mapclimgrid for all region options).The user can also provide a polygon(s) and the mapclimgrid function will clip the grid points outside the specified polygon. Another option is to provide a mask (with dimensions identical to the spatial dimension of the ClimGrid data) which contains NaN and 1.0 and the data inside the ClimGrid struct will be clipped with the mask."
},

{
    "location": "maps/#Timeseries-1",
    "page": "Visualization",
    "title": "Timeseries",
    "category": "section",
    "text": "Plotting timeseries of a given ClimGrid C is simply done by calling plot. This returns the spatial average throughout the time dimension.using ClimateTools\npoly_reg = [[NaN -65 -80 -80 -65 -65];[NaN 42 42 52 52 42]]\n# Extract tasmax variable over specified polygon, between January 1st 1950 and December 31st 2005\nC_hist = load(\"historical.nc\", \"tasmax\", data_units=\"Celsius\", poly=poly_reg, start_date=Date(1950, 01, 01), end_date=Date(2005, 12, 31)))\n# Extract tasmax variable over specified polygon, between January 1st 2006 and December 31st 2090 for emission scenario RCP8.5\nC_future85 = load(\"futureRCP85.nc\", \"tasmax\", data_units=\"Celsius\", poly=poly_reg, start_date=Date(2006, 01, 01), end_date=Date(2090, 12, 31)))\nC = merge(C_hist, C_future)\nind = annualmax(C) # compute annual maximum\nplot(ind)(Image: annualmaxtasmax)Note. Time labels ticks should be improved!The timeserie represent the spatial average of the annual maximum temperature over the following region.mapclimgrid(ind, region = \"QuebecNSP\")(Image: annualmaxtasmax_maps)The map represent the time average over 1950-2090 of the annual maximum temperature."
},

{
    "location": "interface/#",
    "page": "Interface",
    "title": "Interface",
    "category": "page",
    "text": ""
},

{
    "location": "interface/#Interface-1",
    "page": "Interface",
    "title": "Interface",
    "category": "section",
    "text": ""
},

{
    "location": "interface/#Merging-ClimGrid-type-1",
    "page": "Interface",
    "title": "Merging ClimGrid type",
    "category": "section",
    "text": "Sometimes, the timeseries are split among multiple files (e.g. climate models outputs). To obtain the complete timeseries, you can merge 2 ClimGrid. The method is based on AxisArrays merging method and is overloaded for the ClimGrid type.C = merge(C1::ClimGrid, C2::ClimGrid)To merge multiple ClimGrid form an array of files, load has a method that accepts an array of files to merge."
},

{
    "location": "interface/#Operators-1",
    "page": "Interface",
    "title": "Operators",
    "category": "section",
    "text": "Basic statistical functions are overloaded on ClimGrid.mean minimum maximum std varBasic arithmetic operators are also loaded.D = C + 2.0 # will add 2.0 to all elements of C\nD = C::ClimGrid - A::ClimGrid # subtract A from C (useful for climatological difference between a future and historical period\nD = C / A # Ratio of 2 ClimGrids"
},

{
    "location": "interface/#Exporting-1",
    "page": "Interface",
    "title": "Exporting",
    "category": "section",
    "text": "Exporting a ClimGrid to disk to a netCDF format can be done with the write functiond.write(C::ClimGrid, filename::String)"
},

{
    "location": "examples/#",
    "page": "Examples",
    "title": "Examples",
    "category": "page",
    "text": ""
},

{
    "location": "examples/#Examples-1",
    "page": "Examples",
    "title": "Examples",
    "category": "section",
    "text": ""
},

{
    "location": "examples/#Typical-workflow-for-a-climate-scenario-1",
    "page": "Examples",
    "title": "Typical workflow for a climate scenario",
    "category": "section",
    "text": "Note. Climate data can be downloaded at ESGF nodesNote 2. The following example is somewhat convoluted, but it gives an overview of the main steps allowed by ClimateTools package."
},

{
    "location": "examples/#Exploration-1",
    "page": "Examples",
    "title": "Exploration",
    "category": "section",
    "text": "First step before extracting the data is to explore the actual dataset at hand. The function Dataset (reexported from te NCDatasets.jl package) is used to examine the file(s). In this example, the simulation is from the MIROC5 model.Dataset(\"datafile.nc\")\n\nDataset: /path/to/file/tasmax_day_MIROC5_historical_r1i1p1_19900101-19991231.nc\nGroup: /\n\nDimensions\n   time = 3650\n   lat = 128\n   lon = 256\n   bnds = 2\n\nVariables\n  time   (3650)\n    Datatype:    Float64\n    Dimensions:  time\n    Attributes:\n     bounds               = time_bnds\n     units                = days since 1850-1-1\n     calendar             = noleap\n     axis                 = T\n     long_name            = time\n     standard_name        = time\n\n  time_bnds   (2 × 3650)\n    Datatype:    Float64\n    Dimensions:  bnds × time\n\n  lat   (128)\n    Datatype:    Float64\n    Dimensions:  lat\n    Attributes:\n     bounds               = lat_bnds\n     units                = degrees_north\n     axis                 = Y\n     long_name            = latitude\n     standard_name        = latitude\n\n  lat_bnds   (2 × 128)\n    Datatype:    Float64\n    Dimensions:  bnds × lat\n\n  lon   (256)\n    Datatype:    Float64\n    Dimensions:  lon\n    Attributes:\n     bounds               = lon_bnds\n     units                = degrees_east\n     axis                 = X\n     long_name            = longitude\n     standard_name        = longitude\n\n  lon_bnds   (2 × 256)\n    Datatype:    Float64\n    Dimensions:  bnds × lon\n\n  height  \n    Attributes:\n     units                = m\n     axis                 = Z\n     positive             = up\n     long_name            = height\n     standard_name        = height\n\n  tasmax   (256 × 128 × 3650)\n    Datatype:    Float32\n    Dimensions:  lon × lat × time\n    Attributes:\n     standard_name        = air_temperature\n     long_name            = Daily Maximum Near-Surface Air Temperature\n     units                = K\n     original_name        = T2\n     cell_methods         = time: maximum\n     cell_measures        = area: areacella\n     history              = 2011-10-19T12:39:31Z altered by CMOR: Treated scalar dimension: \'height\'. 2011-10-19T12:39:31Z altered by CMOR: replaced missing value flag (-999) with standard missing value (1e+20). 2011-10-19T12:39:31Z altered by CMOR: Inverted axis: lat.\n     coordinates          = height\n     missing_value        = 1.0e20\n     _FillValue           = 1.0e20\n     associated_files     = baseURL: http://cmip-pcmdi.llnl.gov/CMIP5/dataLocation gridspecFile: gridspec_atmos_fx_MIROC5_historical_r0i0p0.nc areacella: areacella_fx_MIROC5_historical_r0i0p0.nc\n\nGlobal attributes\n  institution          = AORI (Atmosphere and Ocean Research Institute, The University of Tokyo, Chiba, Japan), NIES (National Institute for Environmental Studies, Ibaraki, Japan), JAMSTEC (Japan Agency for Marine-Earth Science and Technology, Kanagawa, Japan)\n  institute_id         = MIROC\n  experiment_id        = historical\n  source               = MIROC5 2010 atmosphere: MIROC-AGCM6 (T85L40); ocean: COCO (COCO4.5, 256x224 L50); sea ice: COCO (COCO4.5); land: MATSIRO (MATSIRO, L6); aerosols: SPRINTARS (SPRINTARS 5.00, T85L40)\n  model_id             = MIROC5\n  forcing              = GHG, SA, Oz, LU, Sl, Vl, SS, Ds, BC, MD, OC (GHG includes CO2, N2O, methane, and fluorocarbons; Oz includes OH and H2O2; LU excludes change in lake fraction)\n  parent_experiment_id = piControl\n  parent_experiment_rip = r1i1p1\n  branch_time          = 150015.0\n  contact              = Masahiro Watanabe (hiro@aori.u-tokyo.ac.jp), Seita Emori (emori@nies.go.jp), Masayoshi Ishii (ism@jamstec.go.jp), Masahide Kimoto (kimoto@aori.u-tokyo.ac.jp)\n  references           = Watanabe et al., 2010: Improved climate simulation by MIROC5: Mean states, variability, and climate sensitivity. J. Climate, 23, 6312-6335\n  initialization_method = 1\n  physics_version      = 1\n  tracking_id          = 54e617f1-31a5-47fd-bd57-8736bb7d00ef\n  product              = output\n  experiment           = historical\n  frequency            = day\n  creation_date        = 2011-10-19T12:39:31Z\n  history              = 2011-10-19T12:39:31Z CMOR rewrote data to comply with CF standards and CMIP5 requirements.\n  Conventions          = CF-1.4\n  project_id           = CMIP5\n  table_id             = Table day (26 July 2011) f21c16b785432e6bd3f72e80f2cade49\n  title                = MIROC5 model output prepared for CMIP5 historical\n  parent_experiment    = pre-industrial control\n  modeling_realm       = atmos\n  realization          = 1\n  cmor_version         = 2.7.1\nYou can see the dimensions of the data, as well as the name of the variable(s), in this case \"tasmax\"."
},

{
    "location": "examples/#Extraction-1",
    "page": "Examples",
    "title": "Extraction",
    "category": "section",
    "text": "Now, say you need to create a climate scenario, using a given simulation, over a region defined by the following polygon.poly_reg = [[NaN -65 -80 -80 -65 -65];[NaN 42 42 52 52 42]]\n2×6 Array{Float64,2}:\n NaN  -65.0  -80.0  -80.0  -65.0  -65.0\n NaN   42.0   42.0   52.0   52.0   42.0The extraction of the desired variable can be done with the load function, by providing the polygon.gcmfiles =[\"tasmax_day_MIROC5_historical_r1i1p1_19800101-19891231.nc\",\n\"tasmax_day_MIROC5_historical_r1i1p1_19900101-19991231.nc\",\n\"tasmax_day_MIROC5_historical_r1i1p1_20000101-20091231.nc\"]\n\nmodel = load(gcm_files, \"tasmax\", poly=poly_reg)\nClimGrid struct with data:\n   3-dimensional AxisArray{Float32,3,...} with axes:    \n    :lon, [-78.75, -77.3438, -75.9375, -74.5313, -73.125, -71.7188, -70.3125, -68.9063, -67.5, -66.0938]\n    :lat, [42.7233, 44.1241, 45.5249, 46.9256, 48.3264, 49.7271, 51.1279]\n    :time, Date[1980-01-01, 1980-01-02, 1980-01-03, 1980-01-04, 1980-01-05, 1980-01-06, 1980-01-07, 1980-01-08, 1980-01-09, 1980-01-10  …  2009-12-22, 2009-12-23, 2009-12-24, 2009-12-25, 2009-12-26, 2009-12-27, 2009-12-28, 2009-12-29, 2009-12-30, 2009-12-31]\nAnd data, a 10×7×10950 Array{Float32,3}\nProject: CMIP5\nInstitute: MIROC\nModel: MIROC5\nExperiment: historical\nRun: r1i1p1\nVariable: tasmax\nData units: K\nFrequency: day\nGlobal attributes: Dict{Any,Any} with 27 entries\nFilename: tasmax_day_MIROC5_historical_r1i1p1_19800101-19891231.ncOne possible verification of the extracted data is to map the time-mean data with the mapclimgrid function to see if there is something wrong.mapclimgrid(model, region = \"Quebec\")Which should return the following map.(Image: MIROC5)"
},

{
    "location": "examples/#Sidenote:-merging-files/data-1",
    "page": "Examples",
    "title": "Sidenote: merging files/data",
    "category": "section",
    "text": "Climate data files are usually on the order of multiple GBs and institution generally split a single simulation into multiple files. In order to calculate climatologies, it is thus essential to merge the data into a single structure. The function merge is provided to combine 2 ClimGrid.C = merge(C1, C2) # merge C1 with C2The merge function is useful when you have 2 or 3 files. However, a single simulation can sometimes be splitted into yearly files. Hence, extracting timeseries on climatological timescales can imply loading more than a hundred files just to get a complete timeserie for a given gridpoint. The function load has a method where the 1st positional argument is an Array of strings (as opposed to a single string).C = load(strarray::Array{String,1}, variable::String; poly, start_date::Tuple, end_date::Tuple, data_units::String))This is how the MIROC5 simulation has been loaded."
},

{
    "location": "examples/#Bias-correction-1",
    "page": "Examples",
    "title": "Bias correction",
    "category": "section",
    "text": "An important step in climate scenarios design is to correct the statistical bias of the simulations compared against a chosen reference (more often than not, weather observations). A typical method is to do quantile-quantile mapping between the simulation timeseries and observed timeseries. The function qqmap does so. First step would be to interpolate the simulated field onto the reference grid. Here we use the dataset provided by the Canadian Forest Service (McKenney et al. 2011) for the interpolation step and the bias correction step.McKenney, D. W., Hutchinson, M.F., Papadopol, P., Lawrence, K., Pedlar, J., Campbell, K., Milewska, E., Hopkinson, R., Price, D., Owen, T. (2011). \"Customized spatial climate models for North America.\" Bulletin of American Meteorological Society-BAMS December: 1612-1622.obsfiles = [\"nrcan_canada_daily_tasmax_1950.nc\",\n\"nrcan_canada_daily_tasmax_1951.nc\",\n\"nrcan_canada_daily_tasmax_1952.nc\",\n\"nrcan_canada_daily_tasmax_1953.nc\",\n\"nrcan_canada_daily_tasmax_1954.nc\",\n\"nrcan_canada_daily_tasmax_1955.nc\",\n\"nrcan_canada_daily_tasmax_1956.nc\",\n\"nrcan_canada_daily_tasmax_1957.nc\",\n\"nrcan_canada_daily_tasmax_1958.nc\",\n\"nrcan_canada_daily_tasmax_1959.nc\",\n\"nrcan_canada_daily_tasmax_1960.nc\",\n\"nrcan_canada_daily_tasmax_1961.nc\",\n\"nrcan_canada_daily_tasmax_1962.nc\",\n\"nrcan_canada_daily_tasmax_1963.nc\",\n\"nrcan_canada_daily_tasmax_1964.nc\",\n\"nrcan_canada_daily_tasmax_1965.nc\",\n\"nrcan_canada_daily_tasmax_1966.nc\",\n\"nrcan_canada_daily_tasmax_1967.nc\",\n\"nrcan_canada_daily_tasmax_1968.nc\",\n\"nrcan_canada_daily_tasmax_1969.nc\"]\n\nobs = load(obsfiles, \"tasmax\", poly=poly_reg)mapclimgrid(obs, region = \"Quebec\", titlestr=\"Gridded Obs, 1980-2009\")(Image: NRCAN)"
},

{
    "location": "examples/#Interpolation-/-Regridding-1",
    "page": "Examples",
    "title": "Interpolation / Regridding",
    "category": "section",
    "text": "The interpolation is done with the regrid function. The following command interpolate the values of ClimGrid model onto the grid of ClimGrid obs.modelinterp = regrid(model, obs)\nProgress: 100%|█████████████████████████████████████████| Time: 0:00:38\nClimGrid struct with data:\n   3-dimensional AxisArray{Float64,3,...} with axes:    \n    :lon, Float32[-79.9583, -79.875, -79.7917, -79.7083, -79.625, -79.5417, -79.4583, -79.375, -79.2917, -79.2083  …  -65.7917, -65.7083, -65.625, -65.5417, -65.4583, -65.375, -65.2917, -65.2083, -65.125, -65.0417]\n    :lat, Float32[51.9583, 51.875, 51.7917, 51.7083, 51.625, 51.5417, 51.4583, 51.375, 51.2917, 51.2083  …  42.7917, 42.7083, 42.625, 42.5417, 42.4583, 42.375, 42.2917, 42.2083, 42.125, 42.0417]\n    :time, Date[1980-01-01, 1980-01-02, 1980-01-03, 1980-01-04, 1980-01-05, 1980-01-06, 1980-01-07, 1980-01-08, 1980-01-09, 1980-01-10  …  2009-12-22, 2009-12-23, 2009-12-24, 2009-12-25, 2009-12-26, 2009-12-27, 2009-12-28, 2009-12-29, 2009-12-30, 2009-12-31]\nAnd data, a 180×120×10950 Array{Float64,3}\nProject: CMIP5\nInstitute: MIROC\nModel: MIROC5\nExperiment: historical\nRun: r1i1p1\nVariable: tasmax\nData units: K\nFrequency: day\nGlobal attributes: Dict{Any,Any} with 27 entries\nFilename: tasmax_day_MIROC5_historical_r1i1p1_19800101-19891231.ncjulia> mapclimgrid(modelinterp, region = \"Quebec\", titlestr=\"MIROC5 - Interpolated - 1980-2009\")(Image: MIROC5_INTERP)Notice that there is no new information created here. The interpolation is using Scipy\'s griddata under the hood and is simply a linear interpolation onto the obs grid."
},

{
    "location": "examples/#Quantile-quantile-mapping-1",
    "page": "Examples",
    "title": "Quantile-quantile mapping",
    "category": "section",
    "text": "The high-resolution local information is integrated into ClimGrid modelinterp at the bias correction step. There is a daily transfer function applied on a quantile basis.The call signature is qqmap(obs, ref, fut) where the transfer function is estimated between obs and ref and applied on fut. Note that ref and fut can be the same, as in this example. A typical use-case would be obs and ref covering the same (historical, e.g. 1961-2010) temporal window and fut being a simulation covering a future climatological period (which could be a mix of historic and future, such as 1961-2090). This step is computationally intensive (uses of multiple threads can help here if set by the user).model_qqmap = qqmap(obs, modelinterp, modelinterp)\nProgress: 100%|█████████████████████████████████████████| Time: 0:11:20\nClimGrid struct with data:\n   3-dimensional AxisArray{Float64,3,...} with axes:    \n    :lon, Float32[-79.9583, -79.875, -79.7917, -79.7083, -79.625, -79.5417, -79.4583, -79.375, -79.2917, -79.2083  …  -65.7917, -65.7083, -65.625, -65.5417, -65.4583, -65.375, -65.2917, -65.2083, -65.125, -65.0417]\n    :lat, Float32[51.9583, 51.875, 51.7917, 51.7083, 51.625, 51.5417, 51.4583, 51.375, 51.2917, 51.2083  …  42.7917, 42.7083, 42.625, 42.5417, 42.4583, 42.375, 42.2917, 42.2083, 42.125, 42.0417]\n    :time, Date[1980-01-01, 1980-01-02, 1980-01-03, 1980-01-04, 1980-01-05, 1980-01-06, 1980-01-07, 1980-01-08, 1980-01-09, 1980-01-10  …  2009-12-22, 2009-12-23, 2009-12-24, 2009-12-25, 2009-12-26, 2009-12-27, 2009-12-28, 2009-12-29, 2009-12-30, 2009-12-31]\nAnd data, a 180×120×10950 Array{Float64,3}\nProject: CMIP5\nInstitute: MIROC\nModel: MIROC5\nExperiment: historical\nRun: r1i1p1\nVariable: tasmax\nData units: K\nFrequency: NA\nGlobal attributes: Dict{Any,Any} with 0 entries\nFilename: tasmax_day_MIROC5_historical_r1i1p1_19800101-19891231.ncMapping the results show that the local information is integrated into the model, and that the natural \"mask\" of the observation grid is applied naturally.mapclimgrid(model_qqmap, region = \"Quebec\", titlestr=\"MIROC5 - Interpolated and Quantile-Quantile corrected - 1980-2009\")(Image: MIROC5_QQMAP)Proper assessment of future climate conditions over the specified region would involve replicating these steps for minimally a dozen simulations from multiple models and different emission scenarios (e.g. RCP4.5, RCP8.5, etc.).We can show the effect of bias correction by simply subtracting model_qqmap from modelinterp.mapclimgrid(modelinterp-model_qqmap, region = \"qc\", titlestr=\"MIROC5 - bias correction effect - 1980-2009\", center_cs=true)(Image: MIROC5_effect)"
},

{
    "location": "examples/#Climate-indices-1",
    "page": "Examples",
    "title": "Climate indices",
    "category": "section",
    "text": "Once the climate data is downscaled (interpolation and bias correction) to the proper scale, the user can compute climate indices. For example, annual maximum values of daily maximum temperature could be desired (annualmax).annmax = annualmax(model_qqmap)The return value of climate indices functions are another ClimGrid, but at the yearly scale in the case of annual maximum. Maps and timeseries can be plotted with mapclimgrid and plot respectively.Here\'s the effect of bias correcting on annual maximum values.max_obs = annualmax(obs)\nmax_modelinterp = annualmax(modelinterp)\nmax_modelqqmap = annualmax(model_qqmap)\n\n# Plots\nplot(max_obs, label=\"OBS\")\nplot(max_modelinterp, label=\"MIROC5 - interpolated\")\nplot(max_modelqqmap, label=\"MIROC5 - bias corrected\", titlefig = \"Effect of bias correction on annual maximum values\")(Image: timeseries)"
},

{
    "location": "examples/#Exporting-1",
    "page": "Examples",
    "title": "Exporting",
    "category": "section",
    "text": "Once calculations are done, results can be written to disk to a netCDF file with the write command. Here, we export the annual maximum values of the bias corrected model to the current folder.write(max_modelqqmap, \"annualmax_model_qqmap.nc\")"
},

{
    "location": "functions/#",
    "page": "Index",
    "title": "Index",
    "category": "page",
    "text": ""
},

{
    "location": "functions/#Base.findmax-Tuple{ClimGrid}",
    "page": "Index",
    "title": "Base.findmax",
    "category": "method",
    "text": "findmax(C::ClimGrid; skipnan::Bool=false)\n\nReturn the maximum element of ClimGrid C and its index. If there are multiple maximal elements, then the first one will be returned. If any data element is NaN, this element is returned. The result is in line with max. As climate data can often have NaN values (irregular polygons, missing weather data), the option to skip NaNs is provided as a keyword argument.\n\n\n\n\n\n"
},

{
    "location": "functions/#Base.findmin-Tuple{ClimGrid}",
    "page": "Index",
    "title": "Base.findmin",
    "category": "method",
    "text": "findmin(C::ClimGrid; skipnan::Bool=false)\n\nReturn the minimum element of ClimGrid C and its index. If there are multiple minimal elements, then the first one will be returned. If any data element is NaN, this element is returned. The result is in line with min. As climate data can often have NaN values (irregular polygons, missing weather data), the option to skip NaNs is provided as a keyword argument.\n\n\n\n\n\n"
},

{
    "location": "functions/#Base.maximum-Tuple{ClimGrid}",
    "page": "Index",
    "title": "Base.maximum",
    "category": "method",
    "text": "maximum(A::ClimGrid)\n\nCompute the maximum value of ClimGrid A\n\n\n\n\n\n"
},

{
    "location": "functions/#Base.merge-Tuple{ClimGrid,ClimGrid}",
    "page": "Index",
    "title": "Base.merge",
    "category": "method",
    "text": "merge(A::ClimGrid, B::ClimGrid)\n\nCombines two ClimGrid. Based on the AxisArrays method.\n\n\n\n\n\n"
},

{
    "location": "functions/#Base.minimum-Tuple{ClimGrid}",
    "page": "Index",
    "title": "Base.minimum",
    "category": "method",
    "text": "minimum(A::ClimGrid)\n\nCompute the minimum value of ClimGrid A\n\n\n\n\n\n"
},

{
    "location": "functions/#Base.write-Tuple{ClimGrid,String}",
    "page": "Index",
    "title": "Base.write",
    "category": "method",
    "text": "write(C::ClimGrid, filename::String)\n\nWrite to disk ClimGrid C to netCDF file. Still experimental, some attributes are hardcoded.\n\n\n\n\n\n"
},

{
    "location": "functions/#ClimateTools.annualmax-Tuple{ClimGrid}",
    "page": "Index",
    "title": "ClimateTools.annualmax",
    "category": "method",
    "text": "annualmax(C::ClimGrid)\n\nAnnual maximum of array data.\n\nLet data[i,j] be daily time serie on day i in year j. Extract the highest value for year j.\n\n\n\n\n\n"
},

{
    "location": "functions/#ClimateTools.annualmean-Tuple{ClimGrid}",
    "page": "Index",
    "title": "ClimateTools.annualmean",
    "category": "method",
    "text": "annualmean(C::ClimGrid)\n\nAnnual mean of array data.\n\nLet data[i,j] be daily time serie on day i in year j. Calculate the mean value for year j.\n\n\n\n\n\n"
},

{
    "location": "functions/#ClimateTools.annualmin-Tuple{ClimGrid}",
    "page": "Index",
    "title": "ClimateTools.annualmin",
    "category": "method",
    "text": "annualmin(C::ClimGrid)\n\nAnnual minimum of array data.\n\nLet data[i,j] be daily time serie on day i in year j. Extract the lowest value for year j.\n\n\n\n\n\n"
},

{
    "location": "functions/#ClimateTools.annualsum-Tuple{ClimGrid}",
    "page": "Index",
    "title": "ClimateTools.annualsum",
    "category": "method",
    "text": "annualsum(C::ClimGrid)\n\nAnnual sum of array data.\n\nLet data[i,j] be daily time serie on day i in year j. Sums daily values for year j.\n\n\n\n\n\n"
},

{
    "location": "functions/#ClimateTools.applymask-Tuple{AbstractArray{N,4} where N,AbstractArray{N,2} where N}",
    "page": "Index",
    "title": "ClimateTools.applymask",
    "category": "method",
    "text": "applymask(A::AbstractArray{N, n}, mask::AbstractArray{N, n})\n\nApplies a mask on the array A. Return an AbstractArray{N, n}.\n\n\n\n\n\n"
},

{
    "location": "functions/#ClimateTools.approx_surfacepressure-Tuple{ClimGrid,ClimGrid,ClimGrid}",
    "page": "Index",
    "title": "ClimateTools.approx_surfacepressure",
    "category": "method",
    "text": "  approx_surfacepressure(sealevel_pressure::ClimGrid, orography::ClimGrid, daily_temperature::ClimGrid)\n\nReturns the approximated surface pressure (sp) (Pa) using sea level pressure (psl) (Pa), orography (orog) (m), and daily mean temperature (tas) (K).\n\nsp = psl * 10^x\n\nwhere x = frac-orog18400 * tas  27315\n\n\n\n\n\n"
},

{
    "location": "functions/#ClimateTools.customthresover-Tuple{ClimGrid,Unitful.Quantity}",
    "page": "Index",
    "title": "ClimateTools.customthresover",
    "category": "method",
    "text": "customthresover(C::ClimGrid, thres)\n\ncustomthresover, annual number of days over a specified threshold.\n\nLet TS[i,j] be a daily time serie value on day i in year j. Count the number of days where:\n\nTS[i,j] > thres.\n\nNote. The threshold needs to have units specified. For example:\n\njulia> using Unitful: @u_str, °C julia> thres = 15u\"°C\" 15 °C\n\n\n\n\n\n"
},

{
    "location": "functions/#ClimateTools.customthresunder-Tuple{ClimGrid,Unitful.Quantity}",
    "page": "Index",
    "title": "ClimateTools.customthresunder",
    "category": "method",
    "text": "customthresunder(C::ClimGrid, thres::Quantity)\n\ncustomthresunder, annual number of days under a specified threshold.\n\nLet TS[i,j] be a daily time serie value on day i in year j. Count the number of days where:\n\nTS[i,j] < thres.\n\nNote. The threshold needs to have units specified. For example:\n\njulia> using Unitful: @u_str, °C julia> thres = 15u\"°C\" 15 °C\n\n\n\n\n\n"
},

{
    "location": "functions/#ClimateTools.daymean-Tuple{ClimGrid}",
    "page": "Index",
    "title": "ClimateTools.daymean",
    "category": "method",
    "text": "daymean(C::ClimGrid)\n\nReturns the daily average given a sub-daily ClimGrid.\n\n\n\n\n\n"
},

{
    "location": "functions/#ClimateTools.daysum-Tuple{ClimGrid}",
    "page": "Index",
    "title": "ClimateTools.daysum",
    "category": "method",
    "text": "daysum(C::ClimGrid)\n\nReturns the daily sum given a sub-daily ClimGrid C.\n\n\n\n\n\n"
},

{
    "location": "functions/#ClimateTools.diurnaltemperature-Tuple{ClimGrid,ClimGrid,Float64}",
    "page": "Index",
    "title": "ClimateTools.diurnaltemperature",
    "category": "method",
    "text": "diurnaltemperature(temperatureminimum::ClimGrid, temperaturemaximum::ClimGrid, α::Float64)\n\nReturns an estimation of the diurnal temperature (temperature between 7:00 (7am) and 17:00 (5pm)). The estimation is a linear combination of the daily minimum temperature (temperatureminimum) and daily maximum temperature (temperaturemaximum). The value of α has to be estimated seperatly from observations and depends on the location. The daily max and min must be in the same unit and in Celsius or Kelvin The diurnal temperature returned is in the same units as the daily minimum temperature and daily maximum temperature.\n\nTdiu = α * Tmin + (1 - α) * Tmax\n\n\n\n\n\n"
},

{
    "location": "functions/#ClimateTools.ensemble_max-Tuple{Any}",
    "page": "Index",
    "title": "ClimateTools.ensemble_max",
    "category": "method",
    "text": "ensemble_max(C::ClimGrid...)\n\nReturns the Ensemble maximum of climatological means of ClimGrids C..\n\n\n\n\n\n"
},

{
    "location": "functions/#ClimateTools.ensemble_mean-Tuple{Any}",
    "page": "Index",
    "title": "ClimateTools.ensemble_mean",
    "category": "method",
    "text": "ensemble_mean(C::ClimGrid...)\n\nReturns the Ensemble mean of ClimGrids C..\n\n\n\n\n\n"
},

{
    "location": "functions/#ClimateTools.ensemble_min-Tuple{Any}",
    "page": "Index",
    "title": "ClimateTools.ensemble_min",
    "category": "method",
    "text": "ensemble_min(C::ClimGrid...)\n\nReturns the Ensemble minimum of climatological means of ClimGrids C..\n\n\n\n\n\n"
},

{
    "location": "functions/#ClimateTools.ensemble_std-Tuple{Any}",
    "page": "Index",
    "title": "ClimateTools.ensemble_std",
    "category": "method",
    "text": "ensemble_std(C::ClimGrid...)\n\nReturns the Ensemble standard deviation of climatological means of ClimGrids C..\n\n\n\n\n\n"
},

{
    "location": "functions/#ClimateTools.extractpoly-Tuple{String,Int64}",
    "page": "Index",
    "title": "ClimateTools.extractpoly",
    "category": "method",
    "text": "extractpoly(file::String, n::Int)\n\nReturns the n-th polygon contained in file.\n\n\n\n\n\n"
},

{
    "location": "functions/#ClimateTools.frostdays-Tuple{ClimGrid}",
    "page": "Index",
    "title": "ClimateTools.frostdays",
    "category": "method",
    "text": "frostdays(C::ClimGrid)\n\nFD, Number of frost days: Annual count of days when TN (daily minimum temperature) < 0 Celsius.\n\nLet TN[i,j] be daily minimum temperature on day i in year j. Count the number of days where:\n\nTN[i,j] < 0 Celsius.\n\n\n\n\n\n"
},

{
    "location": "functions/#ClimateTools.get_timevec-Tuple{ClimGrid}",
    "page": "Index",
    "title": "ClimateTools.get_timevec",
    "category": "method",
    "text": "get_timevec(C::ClimGrid)\n\nReturns time vector of ClimGrid C.\n\n\n\n\n\n"
},

{
    "location": "functions/#ClimateTools.icingdays-Tuple{ClimGrid}",
    "page": "Index",
    "title": "ClimateTools.icingdays",
    "category": "method",
    "text": "icingdays(C::ClimGrid)\n\nID, Number of summer days: Annual count of days when TX (daily maximum temperature) < 0 degree Celsius.\n\nLet TX[i,j] be daily maximum temperature on day i in year j. Count the number of days where:\n\nTX[i,j] < 0 Celsius.\n\n\n\n\n\n"
},

{
    "location": "functions/#ClimateTools.inpoly-Tuple{Any,Any}",
    "page": "Index",
    "title": "ClimateTools.inpoly",
    "category": "method",
    "text": "inpoly(p, poly::Matrix)\n\nDetermines if a point is inside a polygon.\n\np – point (x,y) or [x,y]\npoly – polygon vertices [x1 x2 ... xn x1                           y1 y2 ... yn y1]\n\nReturns true if point has an odd winding number.  This should label points as exterior which are inside outcrops.  See test for a test.\n\nAuthor: Github \"Mauro3\" / \"Mauro\"\n\n\n\n\n\n"
},

{
    "location": "functions/#ClimateTools.inpolygrid-Tuple{AbstractArray{N,2} where N,AbstractArray{N,2} where N,AbstractArray{N,2} where N}",
    "page": "Index",
    "title": "ClimateTools.inpolygrid",
    "category": "method",
    "text": "inpolygrid(lon, lat, poly::AbstractArray{N,2} where N)\n\nUsed to test a grid of points. Returns a mask of ones and NaNs of the same size as lon and lat.\n\n\n\n\n\n"
},

{
    "location": "functions/#ClimateTools.load-Tuple{Array{String,1},String}",
    "page": "Index",
    "title": "ClimateTools.load",
    "category": "method",
    "text": "load(files::Array{String,1}, vari::String; poly = ([]), start_date::Date = Date(-4000), end_date::Date = Date(-4000), data_units::String = \"\")\n\nLoads and merge the files contained in the arrar files.\n\n\n\n\n\n"
},

{
    "location": "functions/#ClimateTools.load-Tuple{String,String}",
    "page": "Index",
    "title": "ClimateTools.load",
    "category": "method",
    "text": "load(file::String, vari::String; poly = Array{Float64}([]), start_date::Tuple, end_date::Tuple, data_units::String = \"\")\n\nReturns a ClimGrid type with the data in file of variable vari inside the polygon poly. Metadata is built-in the ClimGrid type, from the netCDF attributes.\n\nInside the ClimgGrid type, the data is stored into an AxisArray data type, with time, longitude/x and latitude/y dimensions.\n\nThe polygon provided should be in the -180, +180 longitude format. If the polygon crosses the International Date Line, the polygon should be splitted in multiple parts (i.e. multi-polygons).\n\nOptions for data_units are for precipitation : \"mm\", which converts the usual \"kg m-2 s-1\" unit found in netCDF files. For temperature : \"Celsius\", which converts the usual \"Kelvin\" unit.\n\nTemporal subsetting can be done by providing start_date and end-date Tuples of length 1 (year), length 3 (year, month, day) or 6 (hour, minute, second).\n\nNote: load uses CF conventions. If you are unable to read the netCDF file with load, the user will need to read it with low-level functions available in NetCDF.jl package or NCDatasets.jl or re-create standartized netCDF files.\n\n\n\n\n\n"
},

{
    "location": "functions/#ClimateTools.load2D-Tuple{String,String}",
    "page": "Index",
    "title": "ClimateTools.load2D",
    "category": "method",
    "text": "load2D(file::String, vari::String; poly=[], data_units::String=\"\")\n\nReturns a 2D array. Should be used for fixed data, such as orography\n\n\n\n\n\n"
},

{
    "location": "functions/#ClimateTools.mapclimgrid-Tuple{ClimGrid}",
    "page": "Index",
    "title": "ClimateTools.mapclimgrid",
    "category": "method",
    "text": "mapclimgrid(C::ClimGrid; region::String=\"auto\", poly, level, mask, caxis, start_date::Tuple, end_date::Tuple, titlestr::String, surface::Symbol, cm::String=\"\", ncolors::Int, center_cs::Bool, filename::String, cs_label::String)\n\nMaps the time-mean average of ClimGrid C. If a filename is provided, the figure is saved in a png format.\n\nOptional keyworkd includes precribed regions (keyword region, see list below), spatial clipping by polygon (keyword poly) or mask (keyword mask, an array of NaNs and 1.0 of the same dimension as the data in ClimGrid C), startdate and enddate. For 4D data, keyword level is used to map a given level (defaults to 1). caxis is used to limit the colorscale. cm is used to manually set the colorscale (see Python documentation for native colorscale keyword), ncolors is used to set the number of color classes (defaults to 12). Set center_cs to true to center the colorscale (useful for divergent results, such as anomalies, positive/negative temprature). cs_label is used for custom colorscale label.\n\nArguments for keyword region (and shortcuts)\n\nEurope (\"EU\")\nNorthAmerica (\"NA\")\nCanada (\"CA\")\nQuebec, QuebecNSP (\"QC\", \"QCNSP\")\nAmericas (\"Ams\")\nWorld, WorldAz, WorldEck4 (\"W\", \"Waz\", \"Weck4\")\nGreenwich (\"Gr\")\n\nArguments for keyword surface\n\n:contour\n:contourf\n:pcolormesh\n\n\n\n\n\n"
},

{
    "location": "functions/#ClimateTools.mapclimgrid-Tuple{}",
    "page": "Index",
    "title": "ClimateTools.mapclimgrid",
    "category": "method",
    "text": "mapclimgrid(; region::String=\"auto\", poly, level, mask, caxis, start_date::Date, end_date::Date)\n\nEmpty map generator, when called without a ClimGrid as the positional argument.\n\n\n\n\n\n"
},

{
    "location": "functions/#ClimateTools.meantemperature-Tuple{ClimGrid,ClimGrid}",
    "page": "Index",
    "title": "ClimateTools.meantemperature",
    "category": "method",
    "text": "meantemperature(temperatureminimum::ClimGrid, temperaturemaximum::ClimGrid)\n\nReturns the daily mean temperature calculated from the maximum and minimum temperature. Daily maximum and minimum temperature must be in the same units. The mean temperature returned is in the same units as the daily minimum temperature and daily maximum temperature.\n\nTmean = fracTmax + Tmin2\n\n\n\n\n\n"
},

{
    "location": "functions/#ClimateTools.meshgrid-Union{Tuple{T}, Tuple{AbstractArray{T,1},AbstractArray{T,1}}} where T",
    "page": "Index",
    "title": "ClimateTools.meshgrid",
    "category": "method",
    "text": "X, Y = meshgrid{T}(vx::AbstractVector{T}, vy::AbstractVector{T})\n\nThis function creates a 2-D mesh-grid in a format consistent with Matlab\'s function meshgrid(). XV and YV are vectors.\n\n\n\n\n\n"
},

{
    "location": "functions/#ClimateTools.monthmean-Tuple{ClimGrid}",
    "page": "Index",
    "title": "ClimateTools.monthmean",
    "category": "method",
    "text": "monthmean(C::ClimGrid)\n\nReturns monthly means of ClimGrid C.\n\n\n\n\n\n"
},

{
    "location": "functions/#ClimateTools.monthsum-Tuple{ClimGrid}",
    "page": "Index",
    "title": "ClimateTools.monthsum",
    "category": "method",
    "text": "monthsum(C::ClimGrid)\n\nReturns monthly sums of ClimGrid C.\n\n\n\n\n\n"
},

{
    "location": "functions/#ClimateTools.ndgrid-Tuple{AbstractArray{T,1} where T}",
    "page": "Index",
    "title": "ClimateTools.ndgrid",
    "category": "method",
    "text": "X, Y = ndgrid(XV, YV)\n\nThis function creates a 2-D mesh-grid in a format consistent with Matlab\'s function ndgrid(). XV and YV are vectors.\n\n\n\n\n\n"
},

{
    "location": "functions/#ClimateTools.periodmean-Tuple{ClimGrid}",
    "page": "Index",
    "title": "ClimateTools.periodmean",
    "category": "method",
    "text": "periodmean(C::ClimGrid; startdate::Tuple, enddate::Tuple)\n\nMean of array data over a given period.\n\n\n\n\n\n"
},

{
    "location": "functions/#ClimateTools.polyfit-Tuple{ClimGrid}",
    "page": "Index",
    "title": "ClimateTools.polyfit",
    "category": "method",
    "text": "polyfit(C::ClimGrid)\n\nReturns an array of the polynomials functions of each grid points contained in ClimGrid C.\n\n\n\n\n\n"
},

{
    "location": "functions/#ClimateTools.polyval-Tuple{ClimGrid,Array{Polynomials.Poly{Float64},2}}",
    "page": "Index",
    "title": "ClimateTools.polyval",
    "category": "method",
    "text": "polyval(C::ClimGrid, polynomial::Array{Poly{Float64}})\n\nReturns a ClimGrid containing the values, as estimated from polynomial function polyn.\n\n\n\n\n\n"
},

{
    "location": "functions/#ClimateTools.prcp1-Tuple{ClimGrid}",
    "page": "Index",
    "title": "ClimateTools.prcp1",
    "category": "method",
    "text": "prcp1(C::ClimGrid)\n\nAnnual number with preciptation >= 1 mm. This function returns a ClimGrid. Input data should be in mm.\n\n\n\n\n\n"
},

{
    "location": "functions/#ClimateTools.qqmap-Tuple{Array{N,1} where N,Array{N,1} where N,Array{N,1} where N,Any,Any,Any,Any}",
    "page": "Index",
    "title": "ClimateTools.qqmap",
    "category": "method",
    "text": "qqmap(obs::Array{N, 1} where N, ref::Array{N, 1} where N, fut::Array{N, 1} where N; method=\"Additive\", detrend=true, window=15, rankn=50, thresnan=0.1, keep_original=false, interp::Function = Linear(), extrap::Function = Flat())\n\nQuantile-Quantile mapping bias correction for single vector. This is a low level function used by qqmap(A::ClimGrid ..), but can work independently.\n\n\n\n\n\n"
},

{
    "location": "functions/#ClimateTools.qqmap-Tuple{ClimGrid,ClimGrid,ClimGrid}",
    "page": "Index",
    "title": "ClimateTools.qqmap",
    "category": "method",
    "text": "qqmap(obs::ClimGrid, ref::ClimGrid, fut::ClimGrid; method=\"Additive\", detrend=true, window::Int=15, rankn::Int=50, thresnan::Float64=0.1, keep_original::Bool=false, interp::Function = Linear(), extrap::Function = Flat())\n\nQuantile-Quantile mapping bias correction. For each julian day of the year (+/- window size), a transfer function is estimated through an empirical quantile-quantile mapping.\n\nThe quantile-quantile transfer function between ref and obs is etimated on a julian day (and grid-point) basis with a moving window around the julian day. Hence, for a given julian day, the transfer function is then applied to the fut dataset for a given julian day.\n\nOptions\n\nmethod::String = \"Additive\" (default) or \"Multiplicative\". Additive is used for most climate variables. Multiplicative is usually bounded variables such as precipitation and humidity.\n\ndetrend::Bool = true (default). A 4th order polynomial is adjusted to the time series and the residuals are corrected with the quantile-quantile mapping.\n\nwindow::Int = 15 (default). The size of the window used to extract the statistical characteristics around a given julian day.\n\nrankn::Int = 50 (default). The number of bins used for the quantile estimations. The quantiles uses by default 50 bins between 0.01 and 0.99. The bahavior between the bins is controlled by the interp keyword argument. The behaviour of the quantile-quantile estimation outside the 0.01 and 0.99 range is controlled by the extrap keyword argument.\n\nthresnan::Float64 = 0.1 (default). The fraction is missing values authorized for the estimation of the quantile-quantile mapping for a given julian days. If there is more than treshnan missing values, the output for this given julian days returns NaNs.\n\nkeep_original::Bool = false (default). If keep_original is set to true, the values are set to the original values in presence of too many NaNs.\n\ninterp = Interpolations.Linear() (default). When the data to be corrected lies between 2 quantile bins, the value of the transfer function is linearly interpolated between the 2 closest quantile estimation. The argument is from Interpolations.jl package.\n\nextrap = Interpolations.Flat() (default). The bahavior of the quantile-quantile transfer function outside the 0.01-0.99 range. Setting it to Flat() ensures that there is no \"inflation problem\" with the bias correction. The argument is from Interpolation.jl package.\n\n\n\n\n\n"
},

{
    "location": "functions/#ClimateTools.regrid-Tuple{ClimGrid,AbstractArray{N,T} where N where T,AbstractArray{N,T} where N where T}",
    "page": "Index",
    "title": "ClimateTools.regrid",
    "category": "method",
    "text": "C = regrid(A::ClimGrid, londest::AbstractArray{N, 1} where N, latdest::AbstractArray{N, 1} where N)A\n\nInterpolate ClimGrid A onto lat-lon grid defined by londest and latdest vector or array. If an array is provided, it is assumed that the grid is curvilinear (not a regular lon-lat grid) and the user needs to provide the dimension vector (\"x\" and \"y\") for such a grid.\n\n\n\n\n\n"
},

{
    "location": "functions/#ClimateTools.regrid-Tuple{ClimGrid,ClimGrid}",
    "page": "Index",
    "title": "ClimateTools.regrid",
    "category": "method",
    "text": "C = regrid(A::ClimGrid, B::ClimGrid; method=\"linear\", min=[], max=[])\n\nInterpolate ClimGrid A onto the lon-lat grid of ClimGrid B, where A and B are ClimGrid. Available methods for interpolation are \"linear\" (default), \"nearest\" and \"cubic\".\n\nMin and max optional keyword are used to constraint the results of the interpolation. For example, interpolating bounded fields can lead to unrealilstic values, such as negative precipitation. In that case, one would use min=0.0 to convert negative precipitation to 0.0.\n\n\n\n\n\n"
},

{
    "location": "functions/#ClimateTools.resample-Tuple{ClimGrid,Int64,Int64}",
    "page": "Index",
    "title": "ClimateTools.resample",
    "category": "method",
    "text": "resample(C::ClimGrid, startmonth::Int64, endmonth::Int64)\n\nReturn a resampled subset of ClimGrid C based on months.\n\n\n\n\n\n"
},

{
    "location": "functions/#ClimateTools.resample-Tuple{ClimGrid,String}",
    "page": "Index",
    "title": "ClimateTools.resample",
    "category": "method",
    "text": "resample(C::ClimGrid, season::String)\n\nReturn a resampled subset of ClimGrid C for a given season. Season options are: \"DJF\" (December-February), \"MAM\" (March-May), \"JJA\" (June-August), \"SON\" (September-November)\n\n\n\n\n\n"
},

{
    "location": "functions/#ClimateTools.shapefile_coords-Tuple{Shapefile.Polygon}",
    "page": "Index",
    "title": "ClimateTools.shapefile_coords",
    "category": "method",
    "text": "shapefile_coords(poly::Shapefile.Polygon)\n\nThis function return the polygons contained in shp.shapes[i]. It returns the x and y coordinates vectors.\n\nSee also shapefile_coords_poly, which returns a polygon that ca be used for data extraction of the load.\n\n\n\n\n\n"
},

{
    "location": "functions/#ClimateTools.shapefile_coords_poly-Tuple{Shapefile.Polygon}",
    "page": "Index",
    "title": "ClimateTools.shapefile_coords_poly",
    "category": "method",
    "text": "shapefile_coords_poly(poly::Shapefile.Polygon)\n\nReturn the polygons contained in shp.shapes[i]. It returns an array containing the polygons.\n\nSee also shapefile_coords, which returns vectors as opposed to array. Returned polygon is consistent with the data extraction of the load function.\n\n\n\n\n\n"
},

{
    "location": "functions/#ClimateTools.spatialsubset-Tuple{ClimGrid,Any}",
    "page": "Index",
    "title": "ClimateTools.spatialsubset",
    "category": "method",
    "text": "spatialsubset(C::ClimGrid, poly::Array{N, 2})\n\nReturns the spatial subset of ClimGrid C. The spatial subset is defined by the polygon poly, defined on a -180, +180 longitude reference.\n\n\n\n\n\n"
},

{
    "location": "functions/#ClimateTools.summerdays-Tuple{ClimGrid}",
    "page": "Index",
    "title": "ClimateTools.summerdays",
    "category": "method",
    "text": "summerdays(C::ClimGrid)\n\nSD, Number of summer days: Annual count of days when TX (daily maximum temperature) > 25 degree Celsius.\n\nLet TX[i,j] be daily maximum temperature on day i in year j. Count the number of days where:\n\nTX[i,j] > 25 Celsius.\n\n\n\n\n\n"
},

{
    "location": "functions/#ClimateTools.temporalsubset-Tuple{ClimGrid,Tuple,Tuple}",
    "page": "Index",
    "title": "ClimateTools.temporalsubset",
    "category": "method",
    "text": "function temporalsubset(C::ClimGrid, startdate::Date, enddate::Date)\n\nReturns the temporal subset of ClimGrid C. The temporal subset is defined by a start and end date.\n\n\n\n\n\n"
},

{
    "location": "functions/#ClimateTools.tropicalnights-Tuple{ClimGrid}",
    "page": "Index",
    "title": "ClimateTools.tropicalnights",
    "category": "method",
    "text": "tropicalnights(C::ClimGrid)\n\nTropicalNights, Number of tropical nights: Annual count of days when TN (daily maximum temperature) > 20 degree Celsius.\n\nLet TN[i,j] be daily minimum temperature on day i in year j. Count the number of days where:\n\nTN[i,j] > 20 Celsius.\n\n\n\n\n\n"
},

{
    "location": "functions/#ClimateTools.vaporpressure-NTuple{4,ClimGrid}",
    "page": "Index",
    "title": "ClimateTools.vaporpressure",
    "category": "method",
    "text": "vaporpressure(specific_humidity::ClimGrid, sealevel_pressure::ClimGrid, orography::ClimGrid, daily_temperature::ClimGrid)\n\nReturns the vapor pressure (vp) (Pa) estimated with the specific humidity (q), the sea level pressure (psl) (Pa), the orography (orog) (m) and the daily mean temperature (tas) (K). An approximation of the surface pressure is first computed by using the sea level pressure, orography and the daily mean temperature (see approx_surfacepressure). Then, vapor pressure is calculated by:\n\nvp = fracq * spq+0622\n\n\n\n\n\n"
},

{
    "location": "functions/#ClimateTools.vaporpressure-Tuple{ClimGrid,ClimGrid}",
    "page": "Index",
    "title": "ClimateTools.vaporpressure",
    "category": "method",
    "text": "vaporpressure(surface_pressure::ClimGrid, specific_humidity::ClimGrid)\n\nReturns the vapor pressure (vp) (Pa) based on the surface pressure (sp) (Pa) and the specific humidity (q).\n\nvp = fracq * spq+0622\n\n\n\n\n\n"
},

{
    "location": "functions/#ClimateTools.wbgt-Tuple{ClimGrid,ClimGrid}",
    "page": "Index",
    "title": "ClimateTools.wbgt",
    "category": "method",
    "text": "wbgt(diurnal_temperature::ClimGrid, vapor_pressure::ClimGrid)\n\nReturns the simplified wet-bulb global temperature (wbgt) (Celsius) calculated using the vapor pressure (Pa) of the day and the estimated mean diurnal temperature (Celsius; temperature between 7:00 (7am) and 17:00 (5pm)).\n\nwbgt = 0567 * Tday + 000393 * vp + 394\n\n\n\n\n\n"
},

{
    "location": "functions/#ClimateTools.yearmonthdayhour-Tuple{AbstractCFDateTime}",
    "page": "Index",
    "title": "ClimateTools.yearmonthdayhour",
    "category": "method",
    "text": "yearmonthdayhour(dt::AbstractCFDateTime) -> (Int64, Int64, Int64, Int64)\n\nSimultaneously return the year, month, day and hour parts of dt.\n\n\n\n\n\n"
},

{
    "location": "functions/#PyPlot.plot-Tuple{ClimGrid}",
    "page": "Index",
    "title": "PyPlot.plot",
    "category": "method",
    "text": "plot(C::ClimGrid, titlefig::String, gridfig::Bool, label::String, color, lw, linestyle)\n\nPlots the spatial average timeserie of ClimGrid C.\n\n\n\n\n\n"
},

{
    "location": "functions/#Statistics.mean-Tuple{ClimGrid}",
    "page": "Index",
    "title": "Statistics.mean",
    "category": "method",
    "text": "mean(A::ClimGrid)\n\nCompute the mean of ClimGrid A\n\n\n\n\n\n"
},

{
    "location": "functions/#Statistics.std-Tuple{ClimGrid}",
    "page": "Index",
    "title": "Statistics.std",
    "category": "method",
    "text": "std(A::ClimGrid)\n\nCompute the standard deviation of ClimGrid A\n\n\n\n\n\n"
},

{
    "location": "functions/#Statistics.var-Tuple{ClimGrid}",
    "page": "Index",
    "title": "Statistics.var",
    "category": "method",
    "text": "var(A::ClimGrid)\n\nCompute the variance of ClimGrid A\n\n\n\n\n\n"
},

{
    "location": "functions/#ClimateTools.ClimGrid",
    "page": "Index",
    "title": "ClimateTools.ClimGrid",
    "category": "type",
    "text": "ClimGrid{A <: AxisArray}\n\nIn-memory representation of Climate Forecast netCDF files.\n\nstruct ClimGrid\n\ndata::AxisArray # Data\n\nlongrid::AbstractArray{N,2} where N # the longitude grid\n\nlatgrid::AbstractArray{N,2} where N # the latitude grid\n\nmsk::Array{N, 2} where N # Data mask (NaNs and 1.0)\n\ngrid_mapping::Dict#{String, Any} # bindings for native grid\n\ndimension_dict::Dict\n\nmodel::String\n\nfrequency::String # Day, month, years\n\nexperiment::String # Historical, RCP4.5, RCP8.5, etc.\n\nrun::String\n\nproject::String # CORDEX, CMIP5, etc.\n\ninstitute::String # UQAM, DMI, etc.\n\nfilename::String # Path of the original file\n\ndataunits::String # Celsius, kelvin, etc.\n\nlatunits::String # latitude coordinate unit\n\nlonunits::String # longitude coordinate unit\n\nvariable::String # Type of variable (i.e. can be the same as \"typeofvar\", but it is changed when calculating indices)\n\ntypeofvar::String # Variable type (e.g. tasmax, tasmin, pr)\n\ntypeofcal::String # Calendar type\n\ntimeattrib::Dict # Time attributes (e.g. days since ... )\n\nvarattribs::Dict # Variable attributes dictionary\n\nglobalattribs::Dict # Global attributes dictionary\n\nend\n\n\n\n\n\n"
},

{
    "location": "functions/#ClimateTools.ClimGrid-Tuple{Any}",
    "page": "Index",
    "title": "ClimateTools.ClimGrid",
    "category": "method",
    "text": "ClimGrid(data; longrid=[], latgrid=[], msk=[], grid_mapping=Dict(), dimension_dict=Dict(), model=\"NA\", frequency=\"NA\", experiment=\"NA\", run=\"NA\", project=\"NA\", institute=\"NA\", filename=\"NA\", dataunits=\"NA\", latunits=\"NA\", lonunits=\"NA\", variable=\"NA\", typeofvar=\"NA\", typeofcal=\"NA\", varattribs=Dict(), globalattribs=Dict())\n\nConstructor of the ClimGrid function. Data is an AxisArray. Everything else is optional, but usually needed for further processing (mapping, interpolation, etc...).\n\nstruct ClimGrid\n\ndata::AxisArray # Data \n\nlongrid::AbstractArray{N,2} where N # the longitude grid \n\nlatgrid::AbstractArray{N,2} where N # the latitude grid \n\nmsk::Array{N, 2} where N # Data mask (NaNs and 1.0) \n\ngrid_mapping::Dict#{String, Any} # bindings for native grid \n\ndimension_dict::Dict\n\nmodel::String\n\nfrequency::String # Day, month, years\n\nexperiment::String # Historical, RCP4.5, RCP8.5, etc.\n\nrun::String\n\nproject::String # CORDEX, CMIP5, etc.\n\ninstitute::String # UQAM, DMI, etc.\n\nfilename::String # Path of the original file\n\ndataunits::String # Celsius, kelvin, etc.\n\nlatunits::String # latitude coordinate unit\n\nlonunits::String # longitude coordinate unit\n\nvariable::String # Type of variable (i.e. can be the same as \"typeofvar\", but it is changed when calculating indices)\n\ntypeofvar::String # Variable type (e.g. tasmax, tasmin, pr)\n\ntypeofcal::String # Calendar type\n\ntimeattrib::Dict # Time attributes (e.g. days since ... )\n\nvarattribs::Dict # Variable attributes dictionary\n\nglobalattribs::Dict # Global attributes dictionary\n\nend\n\n\n\n\n\n"
},

{
    "location": "functions/#ClimateTools.TransferFunction",
    "page": "Index",
    "title": "ClimateTools.TransferFunction",
    "category": "type",
    "text": "TransferFunction(itp::Array, method::String, detrend::Bool)\n\nTransfer function used during quantile-quantile mapping bias correction.\n\n\n\n\n\n"
},

{
    "location": "functions/#Index-1",
    "page": "Index",
    "title": "Index",
    "category": "section",
    "text": "CurrentModule = ClimateToolsModules = [ClimateTools]\nPrivate = false\nOrder   = [:function, :type]"
},

]}
