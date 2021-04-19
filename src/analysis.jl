# """
#     findmax(C::ClimGrid; skipnan::Bool=false)

# Return the maximum element of ClimGrid C and its index. If there are multiple maximal elements, then the first one will be returned. If any data element is NaN, this element is returned. The result is in line with max. As climate data can often have NaN values (irregular polygons, missing weather data), the option to skip NaNs is provided as a keyword argument.

# """
# function Base.findmax(C::ClimGrid; skipnan::Bool=false)
#     # Get data
#     data = C[1].data

#     if skipnan
#         idx = findall(data .== NaNMath.maximum(data))[1]
#         val = data[idx]

#     else
#         idx = argmax(data)
#         val = data[idx]
#     end

#     return val, idx

# end

# """
#     findmin(C::ClimGrid; skipnan::Bool=false)

# Return the minimum element of ClimGrid C and its index. If there are multiple minimal elements, then the first one will be returned. If any data element is NaN, this element is returned. The result is in line with min. As climate data can often have NaN values (irregular polygons, missing weather data), the option to skip NaNs is provided as a keyword argument.

# """
# function Base.findmin(C::ClimGrid; skipnan::Bool=false)
#     # Get data
#     data = C[1].data

#     if skipnan
#         idx = findall(data .== NaNMath.minimum(data))[1]
#         val = data[idx]

#     else
#         idx = argmin(data)
#         val = data[idx]
#     end

#     return val, idx

# end

# """
#     dashboard(C::ClimGrid)
#
# This function returns the summary of ClimGrid C, such as maps of mean, maximum and minimum values in ClimGrid C. The annual cycle, the histogram and estimated probability density function and time series.
# """
# function dashboard(C::ClimGrid; poly = ([]))
#
#     # MAPS, WHILE NICE TO LOOK AT, ARE ALREADY COVERED BY THE FUNCTION "mapclimgrid". Should focus here on annual cycles, and other diagnostics such as annualmax and annualmin
#     # # Mean, maximum and minimum values for each grid points
#     # mean_map = mean(annualmean(C)[1].data, 1)[1, :, :]
#     # FD = AxisArray(mean_map, Axis{:lon}(C[1][Axis{:lon}][:]), Axis{:lat}(C[1][Axis{:lat}][:]))
#     # meanmap = ClimGrid(FD, model = C.model, experiment = C.experiment, run = C.run, filename = C.filename, dataunits = C.dataunits, latunits = C.latunits, lonunits = C.lonunits, variable = string(C.variable, " - temporal mean"), typeofvar = C.typeofvar, typeofcal = C.typeofcal)
#     #
#     # max_map = maximum(annualmax(C)[1].data, 1)[1, :, :]
#     # FD = AxisArray(max_map, Axis{:lon}(C[1][Axis{:lon}][:]), Axis{:lat}(C[1][Axis{:lat}][:]))
#     # maxmap = ClimGrid(FD, model = C.model, experiment = C.experiment, run = C.run, filename = C.filename, dataunits = C.dataunits, latunits = C.latunits, lonunits = C.lonunits, variable = string(C.variable, " - temporal maximum"), typeofvar = C.typeofvar, typeofcal = C.typeofcal)
#     #
#     # min_map = minimum(annualmin(C)[1].data, 1)[1, :, :]
#     # FD = AxisArray(min_map, Axis{:lon}(C[1][Axis{:lon}][:]), Axis{:lat}(C[1][Axis{:lat}][:]))
#     # minmap = ClimGrid(FD, model = C.model, experiment = C.experiment, run = C.run, filename = C.filename, dataunits = C.dataunits, latunits = C.latunits, lonunits = C.lonunits, variable = string(C.variable, " - temporal minimum"), typeofvar = C.typeofvar, typeofcal = C.typeofcal)
#
#     # Spatial-temporal histogram and PDF
#
#     # Annual cycles
#
#
#     # Timeseries (annual mean, max and min)
#     mean_ts = mean(annualmean(C)[1].data,2:3)[:, 1, 1]
#     max_ts = mean(annualmax(C)[1].data,2:3)[:, 1, 1]
#     min_ts = mean(annualmin(C)[1].data,2:3)[:, 1, 1]
#
#
#     # Annual cycle
#
#
#     # Plot everything
#     fig = figure("pyplot_subplot_column",figsize=(10,10))
#     subplot(231) # Create the 1st axis of a 3x1 array of axes
#     title("Mean")
#     subplot(232) # Create the 2nd axis of a 3x1 arrax of axes
#     ax = gca() # Get the handle of the current axis
#     #ax[:set_yscale]("log") # Set the y axis to a logarithmic scale
#     grid("on")
#     ylabel("Log Scale")
#     title("Minimum")
#     subplot(233) # Create the 3rd axis of a 3x1 array of axes
#     ax = gca()
#     ax[:set_xscale]("log") # Set the x axis to a logarithmic scale
#     #xlabel("Log Scale")
#     title("Maximum")
#     fig[:canvas][:draw]() # Update the figure
#     suptitle("3x1 Subplot")
#
#     return true
#
#
# end
