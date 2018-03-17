# TODO tests for annualmean and annualsum

# Prcp1
d = Date(2003,1,1):Date(2005,12,31)

data = vcat(ones(Float64, 365), zeros(Float64, 366), ones(Float64, 365))
Results = Array{Int64, 1}(3); Results[1] = 365; Results[2] = 0; Results[3] = 365;
@test prcp1(data, d) == Results

data= vcat(ones(365,1), zeros(366,1), ones(365))
Results = Array{Int64, 2}(3, 1); Results[1] = 365; Results[2] = 0; Results[3] = 365;
@test prcp1(data, d) == Results

data = Array{Float64,3}(1096, 2, 2)
data[:,1,1] = vcat(ones(365,1), zeros(366,1), ones(365))
data[:,2,1] = vcat(ones(365,1), zeros(366,1), ones(365))
data[:,1,2] = vcat(ones(365,1), zeros(366,1), ones(365))
data[:,2,2] = vcat(ones(365,1), zeros(366,1), ones(365))
Results = Array{Int64, 3}(3, 2, 2); Results[1,1,1] = 365; Results[2,1,1] = 0; Results[3,1,1] = 365; Results[1,2,1] = 365; Results[2,2,1] = 0; Results[3,2,1] = 365; Results[1,1,2] = 365; Results[2,1,2] = 0; Results[3,1,2] = 365; Results[1,2,2] = 365; Results[2,2,2] = 0; Results[3,2,2] = 365;
@test prcp1(data, d) == Results

# ClimGrid based tests
axisdata = AxisArray(data, Axis{:time}(d), Axis{:lon}(1:2), Axis{:lat}(1:2))
C = ClimateTools.ClimGrid(axisdata, variable = "pr")
ind = prcp1(C)
@test ind.data.data == Results

# annualsum
d = Date(2003,1,1):Date(2005,12,31)

data = vcat(ones(Float64, 365), zeros(Float64, 366), ones(Float64, 365))
Results = Array{Int64, 1}(3); Results[1] = 365.; Results[2] = 0.; Results[3] = 365.;
@test annualsum(data, d) == Results

data= vcat(ones(365,1), zeros(366,1), ones(365))
Results = Array{Int64, 2}(3, 1); Results[1] = 365.; Results[2] = 0.; Results[3] = 365.;
@test annualsum(data, d) == Results

data = Array{Float64,3}(1096, 2, 2)
data[:,1,1] = vcat(ones(365,1), zeros(366,1), ones(365))
data[:,2,1] = vcat(ones(365,1), zeros(366,1), ones(365))
data[:,1,2] = vcat(ones(365,1), zeros(366,1), ones(365))
data[:,2,2] = vcat(ones(365,1), zeros(366,1), ones(365))
Results = Array{Int64, 3}(3, 2, 2); Results[1,1,1] = 365.; Results[2,1,1] = 0.; Results[3,1,1] = 365.; Results[1,2,1] = 365.; Results[2,2,1] = 0; Results[3,2,1] = 365.; Results[1,1,2] = 365.; Results[2,1,2] = 0.; Results[3,1,2] = 365.; Results[1,2,2] = 365.; Results[2,2,2] = 0.; Results[3,2,2] = 365.;
@test annualsum(data, d) == Results

# ClimGrid based tests
axisdata = AxisArray(data, Axis{:time}(d), Axis{:lon}(1:2), Axis{:lat}(1:2))
C = ClimateTools.ClimGrid(axisdata, variable = "pr")
ind = annualsum(C)
@test ind.data.data == Results

# annualmean
d = Date(2003,1,1):Date(2005,12,31)

data = vcat(ones(Float64, 365), zeros(Float64, 366), ones(Float64, 365))
Results = Array{Int64, 1}(3); Results[1] = 1.; Results[2] = 0.; Results[3] = 1.;
@test annualmean(data, d) == Results

data= vcat(ones(365,1), zeros(366,1), ones(365))
Results = Array{Int64, 2}(3, 1); Results[1] = 1.; Results[2] = 0.; Results[3] = 1.;
@test annualmean(data, d) == Results

data = Array{Float64,3}(1096, 2, 2)
data[:,1,1] = vcat(ones(365,1), zeros(366,1), ones(365))
data[:,2,1] = vcat(ones(365,1), zeros(366,1), ones(365))
data[:,1,2] = vcat(ones(365,1), zeros(366,1), ones(365))
data[:,2,2] = vcat(ones(365,1), zeros(366,1), ones(365))
Results = Array{Int64, 3}(3, 2, 2); Results[1,1,1] = 1.; Results[2,1,1] = 0.; Results[3,1,1] = 1.; Results[1,2,1] = 1.; Results[2,2,1] = 0; Results[3,2,1] = 1.; Results[1,1,2] = 1.; Results[2,1,2] = 0.; Results[3,1,2] = 1.; Results[1,2,2] = 1.; Results[2,2,2] = 0.; Results[3,2,2] = 1.;
@test annualmean(data, d) == Results

# ClimGrid based tests
axisdata = AxisArray(data, Axis{:time}(d), Axis{:lon}(1:2), Axis{:lat}(1:2))
C = ClimateTools.ClimGrid(axisdata, variable = "pr")
ind = annualmean(C)
@test ind.data.data == Results

# Frostdays
data= vcat(-ones(365), -ones(366),zeros(365))
Results = Array{Int64, 1}(3); Results[1] = 365; Results[2] = 366; Results[3] = 0;
@test frostdays(data, d) == Results
@test icingdays(data, d) == Results


data= vcat(-ones(365,1), -ones(366),zeros(365,1))
Results = Array{Int64, 2}(3, 1); Results[1] = 365; Results[2] = 366; Results[3] = 0;
@test frostdays(data, d) == Results
@test icingdays(data, d) == Results

data = Array{Float64,3}(1096, 2, 2)
data[:,1,1] = vcat(-ones(365,1), -ones(366,1), zeros(365))
data[:,2,1] = vcat(-ones(365,1), -ones(366,1), zeros(365))
data[:,1,2] = vcat(-ones(365,1), -ones(366,1), zeros(365))
data[:,2,2] = vcat(-ones(365,1), -ones(366,1), zeros(365))
Results = Array{Int64, 3}(3, 2, 2); Results[1,1,1] = 365; Results[2,1,1] = 366; Results[3,1,1] = 0; Results[1,2,1] = 365; Results[2,2,1] = 366; Results[3,2,1] = 0; Results[1,1,2] = 365; Results[2,1,2] = 366; Results[3,1,2] = 0; Results[1,2,2] = 365; Results[2,2,2] = 366; Results[3,2,2] = 0;
@test frostdays(data, d) == Results
@test icingdays(data, d) == Results

# ClimGrid based tests
axisdata = AxisArray(data, Axis{:time}(d), Axis{:lon}(1:2), Axis{:lat}(1:2))
C = ClimateTools.ClimGrid(axisdata, variable = "tasmin")
ind = frostdays(C)
@test ind.data.data == Results
C = ClimateTools.ClimGrid(axisdata, variable = "tasmax")
ind = icingdays(C)
@test ind.data.data == Results

# Summer Days
d = Date(2003,1,1):Date(2007,12,31)
data = collect(1.:1826.)
@test summerdays(data, d) == [340, 366, 365, 365, 365]

data = hcat(collect(1.:1826.), collect(1.:1826.))
@test summerdays(data, d) == hcat([340, 366, 365, 365, 365],[340, 366, 365, 365, 365])

data = Array{Float64, 3}(1826, 2, 2)
data[:, 1, 1] = collect(1.:1826.); data[:, 1, 2] = collect(1.:1826.);data[:, 2, 1] = collect(1.:1826.);data[:, 2, 2] = collect(1.:1826.);
Results = Array{Int64, 3}(5, 2, 2); Results[:, 1, 1] = [340, 366, 365, 365, 365]''; Results[:, 1, 2] = [340, 366, 365, 365, 365]''; Results[:, 2, 1] = [340, 366, 365, 365, 365]''; Results[:, 2, 2] = [340, 366, 365, 365, 365]'';
@test summerdays(data, d) == Results

# ClimGrid based tests
axisdata = AxisArray(data, Axis{:time}(d), Axis{:lon}(1:2), Axis{:lat}(1:2))
C = ClimateTools.ClimGrid(axisdata, variable = "tasmax")
ind = summerdays(C)
@test ind.data.data == Results

# Tropical Nights
d = Date(2003,1,1):Date(2007,12,31)
data = collect(1.:1826.)
@test tropicalnights(data, d) == [345, 366, 365, 365, 365]

data = hcat(collect(1.:1826.), collect(1.:1826.))
@test tropicalnights(data, d) == hcat([345, 366, 365, 365, 365],[345, 366, 365, 365, 365])

data = Array{Float64, 3}(1826, 2, 2)
data[:, 1, 1] = collect(1.:1826.); data[:, 1, 2] = collect(1.:1826.);data[:, 2, 1] = collect(1.:1826.);data[:, 2, 2] = collect(1.:1826.);
Results = Array{Int64, 3}(5, 2, 2); Results[:, 1, 1] = [345, 366, 365, 365, 365]''; Results[:, 1, 2] = [345, 366, 365, 365, 365]''; Results[:, 2, 1] = [345, 366, 365, 365, 365]''; Results[:, 2, 2] = [345, 366, 365, 365, 365]'';
@test tropicalnights(data, d) == Results

# ClimGrid based tests
axisdata = AxisArray(data, Axis{:time}(d), Axis{:lon}(1:2), Axis{:lat}(1:2))
C = ClimateTools.ClimGrid(axisdata, variable = "tasmin")
ind = tropicalnights(C)
@test ind.data.data == Results

# Custom thresholds
d = Date(2003,1,1):Date(2007,12,31)
data = collect(1.:1826.)
@test customthresover(data, d, 20) == [345, 366, 365, 365, 365]

data = hcat(collect(1.:1826.), collect(1.:1826.))
@test customthresover(data, d, 20) == hcat([345, 366, 365, 365, 365],[345, 366, 365, 365, 365])

data = Array{Float64, 3}(1826, 2, 2)
data[:, 1, 1] = collect(1.:1826.); data[:, 1, 2] = collect(1.:1826.);data[:, 2, 1] = collect(1.:1826.);data[:, 2, 2] = collect(1.:1826.);
Results = Array{Int64, 3}(5, 2, 2); Results[:, 1, 1] = [345, 366, 365, 365, 365]''; Results[:, 1, 2] = [345, 366, 365, 365, 365]''; Results[:, 2, 1] = [345, 366, 365, 365, 365]''; Results[:, 2, 2] = [345, 366, 365, 365, 365]'';
@test customthresover(data, d, 20) == Results

# ClimGrid based tests
axisdata = AxisArray(data, Axis{:time}(d), Axis{:lon}(1:2), Axis{:lat}(1:2))
C = ClimateTools.ClimGrid(axisdata, variable = "tasmin")
ind = customthresover(C, 20)
@test ind.data.data == Results

d = Date(2003,1,1):Date(2007,12,31)
data = collect(-800.:1025.)
@test customthresunder(data, d, 200) == [365, 366, 269, 0, 0]

data = hcat(collect(-800.:1025.), collect(-800.:1025.))
@test customthresunder(data, d, 200) == hcat([365, 366, 269, 0, 0],[365, 366, 269, 0, 0])

data = Array{Float64, 3}(1826, 2, 2)
data[:, 1, 1] = collect(-800.:1025.); data[:, 1, 2] = collect(-800.:1025.);data[:, 2, 1] = collect(-800.:1025.);data[:, 2, 2] = collect(-800.:1025.);
Results = Array{Int64, 3}(5, 2, 2); Results[:, 1, 1] = [365, 366, 269, 0, 0]''; Results[:, 1, 2] = [365, 366, 269, 0, 0]''; Results[:, 2, 1] = [365, 366, 269, 0, 0]''; Results[:, 2, 2] = [365, 366, 269, 0, 0]'';
@test customthresunder(data, d, 200) == Results

# ClimGrid based tests
axisdata = AxisArray(data, Axis{:time}(d), Axis{:lon}(1:2), Axis{:lat}(1:2))
C = ClimateTools.ClimGrid(axisdata, variable = "tasmin")
ind = customthresunder(C, 200)
@test ind.data.data == Results

# ANNUAL MAXIMUM
d = Date(2003,1,1):Date(2007,12,31)
data = collect(-800.:1025.)
@test annualmax(data, d) == [-436., -70., 295., 660., 1025.]

data = hcat(collect(-800.:1025.), collect(-800.:1025.))
@test annualmax(data, d) == hcat([-436., -70., 295., 660., 1025.],[-436., -70., 295., 660., 1025.])

data = Array{Float64, 3}(1826, 2, 2)
data[:, 1, 1] = collect(-800.:1025.); data[:, 1, 2] = collect(-800.:1025.);data[:, 2, 1] = collect(-800.:1025.);data[:, 2, 2] = collect(-800.:1025.);
Results = Array{Float64, 3}(5, 2, 2); Results[:, 1, 1] = [-436., -70., 295., 660., 1025.]''; Results[:, 1, 2] = [-436., -70., 295., 660., 1025.]''; Results[:, 2, 1] = [-436., -70., 295., 660., 1025.]''; Results[:, 2, 2] = [-436., -70., 295., 660., 1025.]'';
@test annualmax(data, d) == Results



# ClimGrid based tests
axisdata = AxisArray(data, Axis{:time}(d), Axis{:lon}(1:2), Axis{:lat}(1:2))
C = ClimateTools.ClimGrid(axisdata, variable = "tasmin")
ind = annualmax(C)
@test ind.data.data == Results

datetmp = Date(2003,01,01):Date(2008,1,1)
# idx = Dates.monthday(datetmp) .== (2,29)
idx = (Dates.month.(datetmp) .== 2) .& (Dates.day.(datetmp) .== 29)
datetmp = datetmp[.!idx] #creates an Array of dates
data = Array{Float64, 3}(1826, 2, 2)
data[:, 1, 1] = collect(-800.:1025.); data[:, 1, 2] = collect(-800.:1025.);data[:, 2, 1] = collect(-800.:1025.);data[:, 2, 2] = collect(-800.:1025.);
Results = Array{Float64, 3}(6, 2, 2); Results[:, 1, 1] = [-436., -71., 294., 659., 1024., 1025.]''; Results[:, 1, 2] = [-436., -71., 294., 659., 1024., 1025]''; Results[:, 2, 1] = [-436., -71., 294., 659., 1024., 1025.]''; Results[:, 2, 2] = [-436., -71., 294., 659., 1024., 1025.]'';
@test annualmax(data, datetmp) == Results

# ANNUAL MINIMUM
d = Date(2003,1,1):Date(2007,12,31)
data = collect(-800.:1025.)
@test annualmin(data, d) == [-800., -435., -69., 296., 661.]

data = hcat(collect(-800.:1025.), collect(-800.:1025.))
@test annualmin(data, d) == hcat([-800., -435., -69., 296., 661.],[-800., -435., -69., 296., 661.])

data = Array{Float64, 3}(1826, 2, 2)
data[:, 1, 1] = collect(-800.:1025.); data[:, 1, 2] = collect(-800.:1025.);data[:, 2, 1] = collect(-800.:1025.);data[:, 2, 2] = collect(-800.:1025.);
Results = Array{Int64, 3}(5, 2, 2); Results[:, 1, 1] = [-800., -435., -69., 296., 661.]''; Results[:, 1, 2] = [-800., -435., -69., 296., 661.]''; Results[:, 2, 1] = [-800., -435., -69., 296., 661.]''; Results[:, 2, 2] = [-800., -435., -69., 296., 661.]'';
@test annualmin(data, d) == Results

# ClimGrid based tests
axisdata = AxisArray(data, Axis{:time}(d), Axis{:lon}(1:2), Axis{:lat}(1:2))
C = ClimateTools.ClimGrid(axisdata, variable = "tasmin")
ind = annualmin(C)
@test ind.data.data == Results
