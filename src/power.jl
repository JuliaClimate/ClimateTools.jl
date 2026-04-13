"""
    LombScargle.lombscargle(cube::YAXArray, kwargs...)

Compute the Lomb-Scargle periodogram for a given data cube.

# Arguments
- `cube::YAXArray`: The input data cube.
- `kwargs...`: Additional keyword arguments.

# Returns
- The Lomb-Scargle periodogram.

"""
function LombScargle.lombscargle(cube::YAXArray, kwargs...)
    return _xmap_call(
        clombscargle,
        cube;
        reduced_dims=:time,
        output_axes=(Dim{:LombScargle}(["Number of Frequencies", "Period with maximal power", "Maximal Power"]),),
        function_args=(collect(lookup(cube, :time)),),
        outtype=Union{Missing, Float64},
    )
end

"""
    clombscargle(xout, xin, times)

Compute the Lomb-Scargle periodogram for a given time series.

# Arguments
- `xout`: An output array to store the results.
- `xin`: The input time series data.
- `times`: The corresponding time values for the data.

# Details
This function computes the Lomb-Scargle periodogram for a given time series using the `LombScargle` package. It first checks if the input time series has at least 10 non-missing values. If not, it sets the output array `xout` to missing and returns. Otherwise, it calculates the date differences and converts them to integer values. Then, it creates a plan for the Lomb-Scargle algorithm using the date intervals and the time series values. Next, it computes the Lomb-Scargle periodogram and finds the maximum period, maximum power, and the number of significant peaks. Finally, it stores the results in the output array `xout`.

"""
function clombscargle(xout, xin, times)
    ind = .!ismissing.(xin)
    ts = collect(nonmissingtype(eltype(xin)), xin[ind])
    x = times[ind]
    if length(ts) < 10
       xout .= missing
       return
    end
    datediff = Date.(x) .- Date(x[1])
    dateint = getproperty.(datediff, :value)
    pl = LombScargle.plan(dateint, ts)
    pgram = LombScargle.lombscargle(pl)
    lsperiod= findmaxperiod(pgram)
    lspower = findmaxpower(pgram)
    lsnum = Float64(LombScargle.M(pgram))
    perval = isempty(lsperiod) ? missing : Float64(lsperiod[1])
    xout .= [lsnum, perval, Float64(lspower)]
end
