"
  prcp1(data::Array, timevector::StepRange{Date,Base.Dates.Day})

Annual number with preciptation over 1 mm. This function returns a boolean vector. *true* if the data is higher or equal to 1 and *false* otherwise."


function prcp1(data::Array, time::StepRange{Date,Base.Dates.Day})
    # @assert length(size(data)) <= 2 # we want 2 dimensions array
    # return sum(map(x -> !isless(x, 1), data), 1)
    return sum([!isless(istep, 1) for istep in data], 1)
    # out = Array{Bool}(size(data, 1), size(data, 2))

end
