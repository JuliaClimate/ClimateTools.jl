"""
This function returns the summary of ClimGrid C, such as maps of mean, maximum and minimum values in ClimGrid C. The annual cycle, the histogram and estimated probability density function and time series.

dashboard(C::ClimGrid)

"""
function dashboard(C::ClimGrid)

    # Mean, maximum and minimum values for each grid points
    mean_map = mean(annualmean(C)[1].data, 1)[1, :, :]
    FD = AxisArray(mean_map, Axis{:lon}(C[1][Axis{:lon}][:]), Axis{:lat}(C[1][Axis{:lat}][:]))
    meanmap = ClimGrid(FD, model = C.model, experiment = C.experiment, run = C.run, filename = C.filename, dataunits = C.dataunits, latunits = C.latunits, lonunits = C.lonunits, variable = string(C.variable, " - temporal mean"), typeofvar = C.typeofvar, typeofcal = C.typeofcal)

    max_map = maximum(annualmax(C)[1].data, 1)[1, :, :]
    FD = AxisArray(max_map, Axis{:lon}(C[1][Axis{:lon}][:]), Axis{:lat}(C[1][Axis{:lat}][:]))
    maxmap = ClimGrid(FD, model = C.model, experiment = C.experiment, run = C.run, filename = C.filename, dataunits = C.dataunits, latunits = C.latunits, lonunits = C.lonunits, variable = string(C.variable, " - temporal maximum"), typeofvar = C.typeofvar, typeofcal = C.typeofcal)

    min_map = minimum(annualmin(C)[1].data, 1)[1, :, :]
    FD = AxisArray(min_map, Axis{:lon}(C[1][Axis{:lon}][:]), Axis{:lat}(C[1][Axis{:lat}][:]))
    minmap = ClimGrid(FD, model = C.model, experiment = C.experiment, run = C.run, filename = C.filename, dataunits = C.dataunits, latunits = C.latunits, lonunits = C.lonunits, variable = string(C.variable, " - temporal minimum"), typeofvar = C.typeofvar, typeofcal = C.typeofcal)

    # Spatial-temporal histogram and PDF


    # Timeseries (annual mean, max and min)
    mean_ts = mean(annualmean(C)[1].data,2:3)[:, 1, 1]
    max_ts = mean(annualmax(C)[1].data,2:3)[:, 1, 1]
    min_ts = mean(annualmin(C)[1].data,2:3)[:, 1, 1]


    # Annual cycle


    # Plot everything
    fig = figure("pyplot_subplot_column",figsize=(10,10))
    subplot(231) # Create the 1st axis of a 3x1 array of axes
    title("Mean")
    subplot(232) # Create the 2nd axis of a 3x1 arrax of axes
    ax = gca() # Get the handle of the current axis
    #ax[:set_yscale]("log") # Set the y axis to a logarithmic scale
    grid("on")
    ylabel("Log Scale")
    title("Minimum")
    subplot(233) # Create the 3rd axis of a 3x1 array of axes
    ax = gca()
    ax[:set_xscale]("log") # Set the x axis to a logarithmic scale
    #xlabel("Log Scale")
    title("Maximum")
    fig[:canvas][:draw]() # Update the figure
    suptitle("3x1 Subplot")

    return true


end
