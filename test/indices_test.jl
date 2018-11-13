# Prcp1
d = collect(DateTime(2003,1,1):Day(1):DateTime(2005,12,31))

data = Array{Float64,3}(undef, 2, 2,1096)
data[1,1,:] = vcat(ones(365,1), zeros(366,1), ones(365))
data[2,1,:] = vcat(ones(365,1), zeros(366,1), ones(365))
data[1,2,:] = vcat(ones(365,1), zeros(366,1), ones(365))
data[2,2,:] = vcat(ones(365,1), zeros(366,1), ones(365))
Results = Array{Int64, 3}(undef, 2, 2, 3);
Results[1,1,1] = 365;
Results[1,1, 2] = 0;
Results[1,1,3] = 365;
Results[1,2,1] = 365;
Results[2,1,1] = 365;
Results[2,1,2] = 0;
Results[2,1,3] = 365;
Results[1,2,1] = 365;
Results[1,2,2] = 0;
Results[1,2,3] = 365;
Results[2,2,1] = 365;
Results[2,2,2] = 0;
Results[2,2,3] = 365;
# @test prcp1(data, d) == Results

# ClimGrid based tests
axisdata = AxisArray(data, Axis{:lon}(1:2), Axis{:lat}(1:2), Axis{:time}(d))
C = ClimateTools.ClimGrid(axisdata, variable = "pr")
ind = prcp1(C)
@test ind.data.data == Results

# annualsum
data = Array{Float64,3}(undef, 2, 2, 1096)
data[1,1,:] = vcat(ones(365,1), zeros(366,1), ones(365))
data[2,1,:] = vcat(ones(365,1), zeros(366,1), ones(365))
data[1,2,:] = vcat(ones(365,1), zeros(366,1), ones(365))
data[2,2,:] = vcat(ones(365,1), zeros(366,1), ones(365))
Results = Array{Int64, 3}(undef, 2, 2, 3);
Results[1,1,1] = 365.;
Results[1,1,2] = 0.;
Results[1,1,3] = 365.;
Results[2,1,1] = 365.;
Results[2,1,2] = 0;
Results[2,1,3] = 365.;
Results[1,2,1] = 365.;
Results[1,2,2] = 0.;
Results[1,2,3] = 365.;
Results[2,2,1] = 365.;
Results[2,2,2] = 0.;
Results[2,2,3] = 365.;
# @test annualsum(data, d) == Results

# ClimGrid based tests
axisdata = AxisArray(data, Axis{:lon}(1:2), Axis{:lat}(1:2), Axis{:time}(d))
C = ClimateTools.ClimGrid(axisdata, variable = "pr")
ind = annualsum(C)
@test ind.data.data == Results

# annualmean
data = Array{Float64,3}(undef, 2, 2, 1096)
data[1,1,:] = vcat(ones(365,1), zeros(366,1), ones(365))
data[2,1,:] = vcat(ones(365,1), zeros(366,1), ones(365))
data[1,2,:] = vcat(ones(365,1), zeros(366,1), ones(365))
data[2,2,:] = vcat(ones(365,1), zeros(366,1), ones(365))
Results = Array{Int64, 3}(undef, 2, 2, 3);
Results[1,1,1] = 1.;
Results[1,1,2] = 0.;
Results[1,1,3] = 1.;
Results[2,1,1] = 1.;
Results[2,1,2] = 0;
Results[2,1,3] = 1.;
Results[1,2,1] = 1.;
Results[1,2,2] = 0.;
Results[1,2,3] = 1.;
Results[2,2,1] = 1.;
Results[2,2,2] = 0.;
Results[2,2,3] = 1.;
# @test annualmean(data, d) == Results

# ClimGrid based tests
axisdata = AxisArray(data, Axis{:lon}(1:2), Axis{:lat}(1:2), Axis{:time}(d))
C = ClimateTools.ClimGrid(axisdata, variable = "pr")
ind = annualmean(C)
@test ind.data.data == Results

# Frostdays

data = Array{Float64,3}(undef, 2, 2, 1096)
data[1,1,:] = vcat(-ones(365,1), -ones(366,1), zeros(365))
data[2,1,:] = vcat(-ones(365,1), -ones(366,1), zeros(365))
data[1,2,:] = vcat(-ones(365,1), -ones(366,1), zeros(365))
data[2,2,:] = vcat(-ones(365,1), -ones(366,1), zeros(365))
Results = Array{Int64, 3}(undef, 2, 2, 3);
Results[1,1,1] = 365;
Results[1,1,2] = 366;
Results[1,1,3] = 0;
Results[2,1,1] = 365;
Results[2,1,2] = 366;
Results[2,1,3] = 0;
Results[1,2,1] = 365;
Results[1,2,2] = 366;
Results[1,2,3] = 0;
Results[2,2,1] = 365;
Results[2,2,2] = 366;
Results[2,2,3] = 0;
# @test frostdays(data, d) == Results
# @test icingdays(data, d) == Results

# ClimGrid based tests
axisdata = AxisArray(data, Axis{:lon}(1:2), Axis{:lat}(1:2), Axis{:time}(d))
C = ClimateTools.ClimGrid(axisdata, variable = "tasmin")
ind = frostdays(C)
@test ind.data.data == Results
C = ClimateTools.ClimGrid(axisdata, variable = "tasmax")
ind = icingdays(C)
@test ind.data.data == Results

# Summer Days
d = collect(DateTime(2003,1,1):Day(1):DateTime(2007,12,31))

data = Array{Float64, 3}(undef, 2, 2, 1826)
data[1,1,:] = collect(1.0:1826.0); data[1,2,:] = collect(1.0:1826.0); data[2,1,:]=collect(1.0:1826.0); data[2,2,:] = collect(1.0:1826.0);
Results = Array{Int64, 3}(undef, 2, 2, 5);
Results[1, 1, :]=[340, 366, 365, 365, 365]'';
Results[1, 2, :] = [340, 366, 365, 365, 365]'';
Results[2, 1, :] = [340, 366, 365, 365, 365]'';
Results[2, 2, :] = [340, 366, 365, 365, 365]'';
# @test summerdays(data, d) == Results

# ClimGrid based tests
axisdata = AxisArray(data, Axis{:lon}(1:2), Axis{:lat}(1:2), Axis{:time}(d))
C = ClimateTools.ClimGrid(axisdata, variable = "tasmax")
ind = summerdays(C)
@test ind.data.data == Results

data = Array{Float64, 3}(undef, 2, 2, 1826)
data[1,1,:] = collect(1.0:1826.0); data[1,2,:]=collect(1.0:1826.0); data[2,1,:]=collect(1.0:1826.0); data[2,2,:]=collect(1.0:1826.0);
Results = Array{Int64, 3}(undef, 2,2,5);
Results[1,1,:]=[345, 366, 365, 365, 365]'';
Results[1,2,:]=[345, 366, 365, 365, 365]'';
Results[2,1,:]=[345, 366, 365, 365, 365]'';
Results[2,2,:]=[345, 366, 365, 365, 365]'';
# @test tropicalnights(data, d) == Results

# ClimGrid based tests
axisdata = AxisArray(data, Axis{:lon}(1:2), Axis{:lat}(1:2), Axis{:time}(d))
C = ClimateTools.ClimGrid(axisdata, variable = "tasmin")
ind = tropicalnights(C)
@test ind.data.data == Results

# Custom thresholds

data = Array{Float64, 3}(undef, 2,2,1826)
data[1,1,:] = collect(1.0:1826.0); data[1,2,:]=collect(1.0:1826.0); data[2,1,:]=collect(1.0:1826.0); data[2,2,:] = collect(1.0:1826.0);
Results = Array{Int64, 3}(undef, 2,2,5);
Results[1,1,:] = [345, 366, 365, 365, 365]'';
Results[1,2,:] = [345, 366, 365, 365, 365]'';
Results[2,1,:] = [345, 366, 365, 365, 365]'';
Results[2,2,:] = [345, 366, 365, 365, 365]'';
# @test customthresover(data, d, 20) == Results

# ClimGrid based tests
axisdata = AxisArray(data, Axis{:lon}(1:2), Axis{:lat}(1:2), Axis{:time}(d))
C = ClimateTools.ClimGrid(axisdata, variable = "tasmin")
ind = customthresover(C, 20)
@test ind.data.data == Results

data = Array{Float64, 3}(undef, 2,2,1826)
data[1,1,:] = collect(-800.0:1025.0); data[1,2,:] = collect(-800.0:1025.0); data[2,1,:]=collect(-800.0:1025.0); data[2,2,:] = collect(-800.0:1025.0);
Results = Array{Int64, 3}(undef, 2,2,5);
Results[1,1,:]=[365, 366, 269, 0, 0]'';
Results[1,2,:]=[365, 366, 269, 0, 0]'';
Results[2,1,:]=[365, 366, 269, 0, 0]'';
Results[2,2,:]=[365, 366, 269, 0, 0]'';
# @test customthresunder(data, d, 200) == Results

# ClimGrid based tests
axisdata = AxisArray(data, Axis{:lon}(1:2), Axis{:lat}(1:2), Axis{:time}(d))
C = ClimateTools.ClimGrid(axisdata, variable = "tasmin")
ind = customthresunder(C, 200)
@test ind.data.data == Results

# ANNUAL MAXIMUM
data = Array{Float64, 3}(undef, 2,2,1826)
data[1,1,:] = collect(-800.0:1025.0); data[1,2,:] = collect(-800.0:1025.0); data[2,1,:] = collect(-800.0:1025.0); data[2,2,:] = collect(-800.0:1025.0);
Results = Array{Float64, 3}(undef, 2,2,5);
Results[1,1,:]=[-436.0, -70.0, 295.0, 660.0, 1025.0]'';
Results[1,2,:]=[-436.0, -70.0, 295.0, 660.0, 1025.0]'';
Results[2,1,:]=[-436.0, -70.0, 295.0, 660.0, 1025.0]'';
Results[2,2,:]=[-436.0, -70.0, 295.0, 660.0, 1025.0]'';
# @test annualmax(data, d) == Results

# ClimGrid based tests
axisdata = AxisArray(data, Axis{:lon}(1:2), Axis{:lat}(1:2), Axis{:time}(d))
C = ClimateTools.ClimGrid(axisdata, variable = "tasmin")
ind = annualmax(C)
@test ind.data.data == Results

# ANNUAL MINIMUM
data = Array{Float64, 3}(undef, 2,2,1826)
data[1,1,:] = collect(-800.:1025.); data[1,2,:]=collect(-800.:1025.); data[2,1,:]=collect(-800.:1025.); data[2,2,:]=collect(-800.:1025.);
Results = Array{Int64, 3}(undef, 2,2,5);
Results[1,1,:]=[-800.0,-435.0,-69.0,296.0,661.0]'';
 Results[1,2,:]=[-800.0,-435.0,-69.0,296.0,661.0]'';
 Results[2,1,:]=[-800.0,-435.0,-69.0,296.0,661.0]'';
 Results[2,2,:]=[-800.0,-435.0,-69.0,296.0,661.0]'';
# @test annualmin(data, d) == Results

# ClimGrid based tests
axisdata = AxisArray(data, Axis{:lon}(1:2), Axis{:lat}(1:2), Axis{:time}(d))
C = ClimateTools.ClimGrid(axisdata, variable = "tasmin")
ind = annualmin(C)
@test ind.data.data == Results

# Days above 10

data = Array{Float64, 3}(undef, 2, 2, 1826)
data[1,1,:] = collect(1.0:1826.0); data[1,2,:] = collect(1.0:1826.0); data[2,1,:]=collect(1.0:1826.0); data[2,2,:] = collect(1.0:1826.0);
Results = Array{Int64, 3}(undef, 2, 2, 5);
Results[1, 1, :]=[356, 366, 365, 365, 365]'';
Results[1, 2, :] = [356, 366, 365, 365, 365]'';
Results[2, 1, :] = [356, 366, 365, 365, 365]'';
Results[2, 2, :] = [356, 366, 365, 365, 365]'';
# @test summerdays(data, d) == Results

# ClimGrid based tests
axisdata = AxisArray(data, Axis{:lon}(1:2), Axis{:lat}(1:2), Axis{:time}(d))
C = ClimateTools.ClimGrid(axisdata, variable = "tas")
ind = daysabove10(C)
@test ind.data.data == Results

# Period mean
data = Array{Float64, 3}(undef, 2, 2, 1826)
data[1,1,:] = collect(1.0:1826.0); data[1,2,:] = collect(-10.0:1815.0); data[2,1,:] = collect(0.0:1825.0); data[2,2,:] = collect(-1725.0:100);
Results = Array{Float64, 2}(undef, 2, 2);
Results[1, 1]= 746.0;
Results[1, 2] = 735.0;
Results[2, 1] = 745.0;
Results[2, 2] = -980.0;

# ClimGrid based tests
axisdata = AxisArray(data, Axis{:lon}(1:2), Axis{:lat}(1:2), Axis{:time}(d))
C = ClimateTools.ClimGrid(axisdata)
ind = periodmean(C, start_date=(2004, 06, 01), end_date=(2005, 08, 31)) #mean between June 1st 2004 and August 31st 2005
@test ind.data.data == Results


Results = Array{Float64, 2}(undef, 2, 2);
Results[1, 1]= 913.5;
Results[1, 2] = 902.5;
Results[2, 1] = 912.5;
Results[2, 2] = -812.5;
ind = periodmean(C)
@test ind.data.data == Results

# Daymean
raw = collect(0:0.25:1459.75) # 6hr data
d = timedecode(raw, "days since 2000-01-01", "noleap")


data = Array{Float64}(undef, 2, 2, 5840)
data[1, 1, :] = collect(1:5840)
data[1, 2, :] = collect(1:5840)
data[2, 1, :] = collect(1:5840)
data[2, 2, :] = collect(1:5840)

axisdata = AxisArray(data, Axis{:lon}(1:2), Axis{:lat}(1:2), Axis{:time}(d))
C = ClimateTools.ClimGrid(axisdata, frequency="6h")
D = daymean(C)
@test D[1][1,1,1] == 2.5
@test D[1][1,1,end] == 5838.5

D = daysum(C)
@test D[1][1,1,1] == 10.0
@test D[1][1,1,end] == 23354.0
