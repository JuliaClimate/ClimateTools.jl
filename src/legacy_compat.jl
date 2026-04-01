"""
    _cube_with_data(template::YAXArray, data)

Create a YAXArray using the axes of `template` and the provided `data`.
"""
function _cube_with_data(template::YAXArray, data)
    return YAXArray(template.axes, data)
end

function _infer_temperature_unit(data; temperature_unit=:auto)
    if temperature_unit != :auto
        return temperature_unit
    end

    values = _valid_values(data)
    if isempty(values)
        return :K
    end

    return mean(values) > 150 ? :K : :Celsius
end

function _as_kelvin(data; temperature_unit=:auto)
    inferred = _infer_temperature_unit(data; temperature_unit=temperature_unit)
    if inferred in (:Celsius, :C, :degC)
        return data .+ 273.15
    end
    return data
end

function _as_celsius(data; temperature_unit=:auto)
    inferred = _infer_temperature_unit(data; temperature_unit=temperature_unit)
    if inferred in (:K, :Kelvin)
        return data .- 273.15
    end
    return data
end

function _yearly_count(xout, xin; predicate, index_list)
    xout .= NaN
    if !isempty(_valid_values(xin))
        for i in eachindex(index_list)
            values = _valid_values(view(xin, index_list[i]))
            if !isempty(values)
                xout[i] = count(predicate, values)
            end
        end
    end
end

function _yearly_count(cube::YAXArray; predicate)
    time_index = year.(cube.time)
    unique_years = unique(time_index)
    index_list = [findall(==(value), time_index) for value in unique_years]
    indims = InDims("time")
    outdims = OutDims(Dim{:time}(dates_builder_yearmonthday_hardcode(unique_years, imois=7, iday=1)))
    return mapCube(_yearly_count, cube; predicate=predicate, indims=indims, outdims=outdims, index_list=index_list, nthreads=Threads.nthreads())
end

"""
    daymean(cube::YAXArray; shifthour=0)

Daily mean convenience wrapper around `daily_fct`.
"""
function daymean(cube::YAXArray; shifthour=0, kwargs...)
    return daily_fct(cube; fct=mean, shifthour=shifthour, kwargs...)
end

"""
    daysum(cube::YAXArray; shifthour=0)

Daily sum convenience wrapper around `daily_fct`.
"""
function daysum(cube::YAXArray; shifthour=0, kwargs...)
    return daily_fct(cube; fct=sum, shifthour=shifthour, kwargs...)
end

annualmax(cube::YAXArray; kwargs...) = yearly_resample(cube; fct=Base.maximum, kwargs...)
annualmin(cube::YAXArray; kwargs...) = yearly_resample(cube; fct=Base.minimum, kwargs...)
annualmean(cube::YAXArray; kwargs...) = yearly_resample(cube; fct=mean, kwargs...)
annualsum(cube::YAXArray; kwargs...) = yearly_resample(cube; fct=Base.sum, kwargs...)

prcp1(cube::YAXArray; threshold=1) = _yearly_count(cube; predicate=x -> x >= threshold)
frostdays(cube::YAXArray) = _yearly_count(cube; predicate=x -> x < 0)
summerdays(cube::YAXArray; threshold=25) = _yearly_count(cube; predicate=x -> x > threshold)
icingdays(cube::YAXArray) = _yearly_count(cube; predicate=x -> x < 0)
tropicalnights(cube::YAXArray; threshold=20) = _yearly_count(cube; predicate=x -> x > threshold)
customthresover(cube::YAXArray, threshold) = _yearly_count(cube; predicate=x -> x > threshold)
customthresunder(cube::YAXArray, threshold) = _yearly_count(cube; predicate=x -> x < threshold)

"""
    meantemperature(temperatureminimum::YAXArray, temperaturemaximum::YAXArray)

Return the mean daily temperature from daily min and max temperature cubes.
"""
function meantemperature(temperatureminimum::YAXArray, temperaturemaximum::YAXArray)
    data = (Array(temperatureminimum) .+ Array(temperaturemaximum)) ./ 2
    return _cube_with_data(temperatureminimum, data)
end

"""
    diurnaltemperature(temperatureminimum::YAXArray, temperaturemaximum::YAXArray, alpha::Real)

Estimate daytime temperature from daily min and max temperature cubes.
"""
function diurnaltemperature(temperatureminimum::YAXArray, temperaturemaximum::YAXArray, alpha::Real)
    data = (alpha .* Array(temperatureminimum)) .+ ((1 - alpha) .* Array(temperaturemaximum))
    return _cube_with_data(temperatureminimum, data)
end

"""
    approx_surfacepressure(sealevel_pressure::YAXArray, orography::YAXArray, daily_temperature::YAXArray; temperature_unit=:auto)

Approximate surface pressure from sea-level pressure, orography and temperature.
"""
function approx_surfacepressure(sealevel_pressure::YAXArray, orography::YAXArray, daily_temperature::YAXArray; temperature_unit=:auto)
    temperature_kelvin = _as_kelvin(Array(daily_temperature); temperature_unit=temperature_unit)
    exponent = (-1.0 .* Array(orography)) ./ (18400.0 .* temperature_kelvin ./ 273.15)
    data = Array(sealevel_pressure) .* (10.0 .^ exponent)
    return _cube_with_data(sealevel_pressure, data)
end

"""
    vaporpressure(specific_humidity::YAXArray, surface_pressure::YAXArray)

Return vapor pressure from specific humidity and surface pressure.
"""
function vaporpressure(specific_humidity::YAXArray, surface_pressure::YAXArray)
    data = (Array(specific_humidity) .* Array(surface_pressure)) ./ (Array(specific_humidity) .+ 0.622)
    return _cube_with_data(surface_pressure, data)
end

"""
    vaporpressure(specific_humidity::YAXArray, sealevel_pressure::YAXArray, orography::YAXArray, daily_temperature::YAXArray; temperature_unit=:auto)

Approximate surface pressure and derive vapor pressure.
"""
function vaporpressure(specific_humidity::YAXArray, sealevel_pressure::YAXArray, orography::YAXArray, daily_temperature::YAXArray; temperature_unit=:auto)
    surface_pressure = approx_surfacepressure(sealevel_pressure, orography, daily_temperature; temperature_unit=temperature_unit)
    return vaporpressure(specific_humidity, surface_pressure)
end

"""
    wbgt(mean_temperature::YAXArray, vapor_pressure::YAXArray; temperature_unit=:auto)

Return the simplified wet-bulb globe temperature in Celsius.
"""
function wbgt(mean_temperature::YAXArray, vapor_pressure::YAXArray; temperature_unit=:auto)
    temperature_celsius = _as_celsius(Array(mean_temperature); temperature_unit=temperature_unit)
    data = (0.567 .* temperature_celsius) .+ (0.00393 .* Array(vapor_pressure)) .+ 3.94
    return _cube_with_data(mean_temperature, data)
end

function timeresolution(timevec::AbstractVector)
    if length(timevec) <= 1
        return "N/A"
    end

    timediff = Base.diff(timevec)[1]
    if timediff == 1.0 || timediff == 1
        return "24h"
    elseif round(timediff, digits=5) == round(12 / 24, digits=5)
        return "12h"
    elseif round(timediff, digits=5) == round(6 / 24, digits=5)
        return "6h"
    elseif round(timediff, digits=5) == round(3 / 24, digits=5)
        return "3h"
    elseif round(timediff, digits=5) == round(1 / 24, digits=5)
        return "1h"
    elseif round(timediff, digits=5) in (365.0, 366.0, 360.0)
        return "Yearly"
    elseif timediff in (28.0, 29.0, 30.0, 31.0)
        return "Monthly"
    end

    return "N/A"
end

function pr_timefactor(rez::String)
    if rez == "24h"
        return 86400.0
    elseif rez == "12h"
        return 43200.0
    elseif rez == "6h"
        return 21600.0
    elseif rez == "3h"
        return 10800.0
    elseif rez == "1h"
        return 3600.0
    end
    return 1.0
end

function daymean_factor(rez::String)
    if rez == "24h" || rez == "day" || rez == "daily"
        return 1
    elseif rez == "12h"
        return 2
    elseif rez == "6h"
        return 4
    elseif rez == "3h"
        return 8
    elseif rez == "1h"
        return 24
    end
    return 1
end