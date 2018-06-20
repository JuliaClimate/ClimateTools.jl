replstr(x) = sprint((io,x) -> show(IOContext(io, :limit => true, :displaysize => (24, 80)), MIME("text/plain"), x), x)
showstr(x) = sprint((io,x) -> show(IOContext(io, :limit => true, :displaysize => (24, 80)), x), x)


# test that load return a ClimGrid type
file1 = joinpath(dirname(@__FILE__), "data", "sresa1b_ncar_ccsm3-example.nc")
file2 = joinpath(dirname(@__FILE__), "data", "sresa1b_ncar_ccsm3-example.nc")
files = [file1, file2]
C = load(files, "tas")
filenc = joinpath(dirname(@__FILE__), "data", "sresa1b_ncar_ccsm3-example.nc")
C = load(filenc, "tas")

@test replstr(C) == "ClimGrid struct with data:\n   3-dimensional AxisArray{Float32,3,...} with axes:\n    :lon, Float32[-180.0, -178.594, -177.188, -175.781, -174.375, -172.969, -171.563, -170.156, -168.75, -167.344  …  165.938, 167.344, 168.75, 170.156, 171.563, 172.969, 174.375, 175.781, 177.188, 178.594]\n    :lat, Float32[-88.9277, -87.5387, -86.1415, -84.7424, -83.3426, -81.9425, -80.5421, -79.1417, -77.7412, -76.3406  …  76.3406, 77.7412, 79.1417, 80.5421, 81.9425, 83.3426, 84.7424, 86.1415, 87.5387, 88.9277]\n    :time, DateTime[2000-05-15T00:00:00]\nAnd data, a 256×128×1 Array{Float32,3}\nProject: IPCC Fourth Assessment\nInstitute: N/A\nModel: N/A\nExperiment: 720 ppm stabilization experiment (SRESA1B)\nRun: N/A\nVariable: tas\nData units: K\nFrequency: N/A\nGlobal attributes: Dict{Any,Any} with 18 entries\nFilename: /home/proy/.julia/v0.6/ClimateTools/test/data/sresa1b_ncar_ccsm3-example.nc"

@test showstr(C) == "ClimateTools.ClimGrid{AxisArrays.AxisArray{Float32,3,Array{Float32,3},Tuple{AxisArrays.Axis{:lon,Array{Float32,1}},AxisArrays.Axis{:lat,Array{Float32,1}},AxisArrays.Axis{:time,Array{DateTime,1}}}}}(Float32[219.223 226.089 … 265.904 265.964; 219.247 226.28 … 265.876 265.948; … ; 219.05 225.677 … 265.957 265.994; 219.214 225.876 … 265.93 265.979], Float32[-180.0 -180.0 … -180.0 -180.0; -178.594 -178.594 … -178.594 -178.594; … ; 177.188 177.188 … 177.188 177.188; 178.594 178.594 … 178.594 178.594], Float32[-88.9277 -87.5387 … 87.5387 88.9277; -88.9277 -87.5387 … 87.5387 88.9277; … ; -88.9277 -87.5387 … 87.5387 88.9277; -88.9277 -87.5387 … 87.5387 88.9277], [1.0 1.0 … 1.0 1.0; 1.0 1.0 … 1.0 1.0; … ; 1.0 1.0 … 1.0 1.0; 1.0 1.0 … 1.0 1.0], Dict(\"grid_mapping_name\"=>\"Regular_longitude_latitude\"), Dict(\"lat\"=>\"lat\",\"lon\"=>\"lon\"), \"N/A\", \"N/A\", \"720 ppm stabilization experiment (SRESA1B)\", \"N/A\", \"IPCC Fourth Assessment\", \"N/A\", \"/home/proy/.julia/v0.6/ClimateTools/test/data/sresa1b_ncar_ccsm3-example.nc\", \"K\", \"degrees_north\", \"degrees_east\", \"tas\", \"tas\", \"noleap\", Dict{String,Any}(Pair{String,Any}(\"cell_methods\", \"time: mean (interval: 1 month)\"),Pair{String,Any}(\"long_name\", \"air_temperature\"),Pair{String,Any}(\"original_units\", \"K\"),Pair{String,Any}(\"history\", \"Added height coordinate\"),Pair{String,Any}(\"standard_name\", \"air_temperature\"),Pair{String,Any}(\"comment\", \"Created using NCL code CCSM_atmm_2cf.ncl on\\n machine eagle163s\"),Pair{String,Any}(\"original_name\", \"TREFHT\"),Pair{String,Any}(\"coordinates\", \"height\"),Pair{String,Any}(\"_FillValue\", 1.0f20),Pair{String,Any}(\"missing_value\", 1.0f20)…), Dict{Any,Any}(Pair{Any,Any}(\"CVS_Id\", \"\\\$Id\\\$\"),Pair{Any,Any}(\"model_name_english\", \"NCAR CCSM\"),Pair{Any,Any}(\"creation_date\", \"\"),Pair{Any,Any}(\"acknowledgment\", \" Any use of CCSM data should acknowledge the contribution\\n of the CCSM project and CCSM sponsor agencies with the \\n following citation:\\n 'This research uses data provided by the Community Climate\\n System Model project (www.ccsm.ucar.edu), supported by the\\n Directorate for Geosciences of the National Science Foundation\\n and the Office of Biological and Environmental Research of\\n the U.S. Department of Energy.'\\nIn addition, the words 'Community Climate System Model' and\\n 'CCSM' should be included as metadata for webpages referencing\\n work using CCSM data or as keywords provided to journal or book\\npublishers of your manuscripts.\\nUsers of CCSM data accept the responsibility of emailing\\n citations of publications of research using CCSM data to\\n ccsm@ucar.edu.\\nAny redistribution of CCSM data must include this data\\n acknowledgement statement.\"),Pair{Any,Any}(\"history\", \"Tue Oct 25 15:08:51 2005: ncks -O -x -v va -m sresa1b_ncar_ccsm3_0_run1_200001.nc sresa1b_ncar_ccsm3_0_run1_200001.nc\\nTue Oct 25 15:07:21 2005: ncks -d time,0 sresa1b_ncar_ccsm3_0_run1_200001_201912.nc sresa1b_ncar_ccsm3_0_run1_200001.nc\\nTue Oct 25 13:29:43 2005: ncks -d time,0,239 sresa1b_ncar_ccsm3_0_run1_200001_209912.nc /var/www/html/tmp/sresa1b_ncar_ccsm3_0_run1_200001_201912.nc\\nThu Oct 20 10:47:50 2005: ncks -A -v va /data/brownmc/sresa1b/atm/mo/va/ncar_ccsm3_0/run1/sresa1b_ncar_ccsm3_0_run1_va_200001_209912.nc /data/brownmc/sresa1b/atm/mo/tas/ncar_ccsm3_0/run1/sresa1b_ncar_ccsm3_0_run1_200001_209912.nc\\nWed Oct 19 14:55:04 2005: ncks -F -d time,01,1200 /data/brownmc/sresa1b/atm/mo/va/ncar_ccsm3_0/run1/sresa1b_ncar_ccsm3_0_run1_va_200001_209912.nc /data/brownmc/sresa1b/atm/mo/va/ncar_ccsm3_0/run1/sresa1b_ncar_ccsm3_0_run1_va_200001_209912.nc\\nWed Oct 19 14:53:28 2005: ncrcat /data/brownmc/sresa1b/atm/mo/va/ncar_ccsm3_0/run1/foo_05_1200.nc /data/brownmc/sresa1b/atm/mo/va/ncar_ccsm3_0/run1/foo_1192_1196.nc /data/brownmc/sresa1b/atm/mo/va/ncar_ccsm3_0/run1/sresa1b_ncar_ccsm3_0_run1_va_200001_209912.nc\\nWed Oct 19 14:50:38 2005: ncks -F -d time,05,1200 /data/brownmc/sresa1b/atm/mo/va/ncar_ccsm3_0/run1/va_A1.SRESA1B_1.CCSM.atmm.2000-01_cat_2099-12.nc /data/brownmc/sresa1b/atm/mo/va/ncar_ccsm3_0/run1/foo_05_1200.nc\\nWed Oct 19 14:49:45 2005: ncrcat /data/brownmc/sresa1b/atm/mo/va/ncar_ccsm3_0/run1/va_A1.SRESA1B_1.CCSM.atmm.2000-01_cat_2079-12.nc /data/brownmc/sresa1b/atm/mo/va/ncar_ccsm3_0/run1/va_A1.SRESA1B_1.CCSM.atmm.2080-01_cat_2099-12.nc /data/brownmc/sresa1b/atm/mo/va/ncar_ccsm3_0/run1/va_A1.SRESA1B_1.CCSM.atmm.2000-01_cat_2099-12.nc\\nCreated from CCSM3 case b30.040a\\n by wgstrand@ucar.edu\\n on Wed Nov 17 14:12:57 EST 2004\\n \\n For all data, added IPCC requested metadata\"),Pair{Any,Any}(\"table_id\", \"Table A1\"),Pair{Any,Any}(\"prg_ID\", \"Source file unknown Version unknown Date unknown\"),Pair{Any,Any}(\"project_id\", \"IPCC Fourth Assessment\"),Pair{Any,Any}(\"comment\", \"This simulation was initiated from year 2000 of \\n CCSM3 model run b30.030a and executed on \\n hardware cheetah.ccs.ornl.gov. The input external forcings are\\nozone forcing    : A1B.ozone.128x64_L18_1991-2100_c040528.nc\\naerosol optics   : AerosolOptics_c040105.nc\\naerosol MMR      : AerosolMass_V_128x256_clim_c031022.nc\\ncarbon scaling   : carbonscaling_A1B_1990-2100_c040609.nc\\nsolar forcing    : Fixed at 1366.5 W m-2\\nGHGs             : ghg_ipcc_A1B_1870-2100_c040521.nc\\nGHG loss rates   : noaamisc.r8.nc\\nvolcanic forcing : none\\nDMS emissions    : DMS_emissions_128x256_clim_c040122.nc\\noxidants         : oxid_128x256_L26_clim_c040112.nc\\nSOx emissions    : SOx_emissions_A1B_128x256_L2_1990-2100_c040608.nc\\n Physical constants used for derived data:\\n Lv (latent heat of evaporation): 2.501e6 J kg-1\\n Lf (latent heat of fusion     ): 3.337e5 J kg-1\\n r[h2o] (density of water      ): 1000 kg m-3\\n g2kg   (grams to kilograms    ): 1000 g kg-1\\n \\n Integrations were performed by NCAR and CRIEPI with support\\n and facilities provided by NSF, DOE, MEXT and ESC/JAMSTEC.\"),Pair{Any,Any}(\"cmd_ln\", \"bds -x 256 -y 128 -m 23 -o /data/zender/data/dst_T85.nc\")…))"

@test load(filenc, "tas", data_units = "Celsius")[2] == "Celsius"
@test load(filenc, "pr", data_units = "mm")[2] == "mm"
@test typeof(load(filenc, "tas")) == ClimateTools.ClimGrid{AxisArrays.AxisArray{Float32,3,Array{Float32,3},Tuple{AxisArrays.Axis{:lon,Array{Float32,1}},AxisArrays.Axis{:lat,Array{Float32,1}},AxisArrays.Axis{:time,Array{DateTime,1}}}}}

@test typeof(ClimateTools.buildtimevec(filenc, "24h")) == Array{DateTime, 1}

fileorog = joinpath(dirname(@__FILE__), "data", "orog_fx_GFDL-ESM2G_historicalMisc_r0i0p0.nc")
orog = load2D(fileorog, "orog")
@test size(orog[1]) == (144, 90)
@test orog.typeofvar == "orog"
status, figh = mapclimgrid(orog); @test status == true; PyPlot.close();

filename = joinpath(dirname(@__FILE__), "data", "SudQC_GCM.shp")
polyshp = read(filename,Shapefile.Handle)
P = shapefile_coords_poly(polyshp.shapes[1])
orog = load2D(fileorog, "orog", poly = P)
@test size(orog[1]) == (13, 8)
@test orog.typeofvar == "orog"

# Time units
units = NetCDF.ncgetatt(filenc, "time", "units") # get starting date
m = match(r"(\d+)[-.\/](\d+)[-.\/](\d+)", units, 1) # match a date from string
daysfrom = m.match # get only the date ()"yyyy-mm-dd" format)
initDate = DateTime(daysfrom, "yyyy-mm-dd")
timeRaw = floor.(NetCDF.ncread(filenc, "time"))
@test ClimateTools.sumleapyear(initDate::DateTime, timeRaw) == 485

# INTERFACE
# B = vcat(C, C)
# @test size(B.data) == (256, 128, 2) # vcat does not look at dimensions
B = merge(C, C)
@test size(B.data) == (256, 128, 1) # C being similar, they should not add up, as opposed to vcat
# Operators +, -, *, /
B = C + C; @test B[1].data[1, 1, 1] == 438.4457f0
B = C * C; @test B[1].data[1, 1, 1] == 48058.66f0
B = C / C; @test B[1].data[1, 1, 1] == 1.0f0
B = C - C; @test B[1].data[1, 1, 1] == 0.0f0
B = C - 1.0; @test B[1].data[1, 1, 1] == 218.2228546142578
B = C - 1; @test B[1].data[1, 1, 1] == 218.22285f0
B = C / 2; @test B[1].data[1, 1, 1] == 109.61143f0
B = C / 2.2; @test B[1].data[1, 1, 1] == 99.6467520973899
B = C * 2; @test B[1].data[1, 1, 1] == 438.4457f0
B = C * 2.2; @test B[1].data[1, 1, 1] == 482.2902801513672

@test mean(C) == 278.6421f0
@test maximum(C) == 309.09613f0
@test minimum(C) == 205.24321f0
@test std(C) == 21.92836f0
@test round(var(C), 3) == 480.853f0

# @test typeof(show(C)) == Dict{Any, Any}
@test typeof(C[1].data) == Array{Float64,3} || typeof(C[1].data) == Array{Float32,3}
@test C[2] == "K"
@test C[3] == "N/A"
@test C[4] == "720 ppm stabilization experiment (SRESA1B)"
@test C[5] == "N/A"
@test C[6] == "degrees_east"
@test C[7] == "degrees_north"
@test C[8] == joinpath(dirname(@__FILE__), "data", "sresa1b_ncar_ccsm3-example.nc")
@test C[9] == "tas"
@test annualmax(C)[9] == "annualmax"
@test C[10] == "tas"
@test C[11] == "noleap"
@test typeof(C[12]) == Dict{Any, Any}
@test C[12]["project_id"] == "IPCC Fourth Assessment"
@test_throws ErrorException C[13]
@test annualmax(C)[10] == "tas"
@test size(C) == (21, )
@test size(C, 1) == 21
@test length(C) == 21
@test endof(C) == 21
@test_throws ErrorException C[end]
@test ndims(C) == 1


# Spatial subset
filename = joinpath(dirname(@__FILE__), "data", "SudQC_GCM.shp")
filenc = joinpath(dirname(@__FILE__), "data", "sresa1b_ncar_ccsm3-example.nc")
polyshp = read(filename,Shapefile.Handle)
x, y = shapefile_coords(polyshp.shapes[1])
P = [x y]
P = P'
C = load(filenc, "tas")
Csub = spatialsubset(C, P)
@test size(Csub[1]) == (23, 12, 1)
@test Csub[1][1, 1, 1] == 294.6609f0
Csub = spatialsubset(C, P')
@test size(Csub[1]) == (23, 12, 1)
@test Csub[1][1, 1, 1] == 294.6609f0
C = load(filenc, "ua")
Csub = spatialsubset(C, P)
@test size(Csub[1]) == (23, 12, 17, 1)
@test Csub[1][12, 1, 1, 1] == 6.658482f0
@test isnan(Csub[1][1, 1, 1, 1])

poly= [[NaN 10 -10 -10 10 10];[NaN -10 -20 10 10 -10]] # meridian test
C = load(filenc, "tas", poly=poly)

# Spatial subset
C = load(filenc, "tas")
Csub = temporalsubset(C, (2000, 05, 15), (2000, 05, 15))
@test Csub[1][1, 1, 1] == 219.22285f0
@test Csub[1][Axis{:time}][1] == DateTime(2000, 05, 15)
B = load(filenc, "tas", start_date=(2000, 05, 15), end_date=(2000, 05, 15))
@test B[1] == C[1]

# Time resolution
timevec = [1, 2, 3]
@test ClimateTools.timeresolution(timevec) == "24h"
timevec = [1.0, 1.5, 2.0]
@test ClimateTools.timeresolution(timevec) == "12h"
timevec = [1.25, 1.5, 1.75]
@test ClimateTools.timeresolution(timevec) == "6h"
timevec = [1.125, 1.25, 1.375]
@test ClimateTools.timeresolution(timevec) == "3h"
timevec = NetCDF.ncread(filenc, "time")
@test ClimateTools.timeresolution(timevec) == "N/A"



# MESHGRID
YV = [1, 2, 3]
XV = [4, 5, 6]
@test meshgrid(XV, YV) == ([4 5 6; 4 5 6; 4 5 6], [1 1 1; 2 2 2; 3 3 3])
@test meshgrid(XV) == ([4 5 6; 4 5 6; 4 5 6], [4 4 4; 5 5 5; 6 6 6])
Q = Array{Int64, 3}(3, 3, 3)
R = Array{Int64, 3}(3, 3, 3)
S = Array{Int64, 3}(3, 3, 3)
Q[:, :, 1] = [4 5 6; 4 5 6; 4 5 6]
Q[:, :, 2] = [4 5 6; 4 5 6; 4 5 6]
Q[:, :, 3] = [4 5 6; 4 5 6; 4 5 6]
R[:, :, 1] = [1 1 1; 2 2 2; 3 3 3]
R[:, :, 2] = [1 1 1; 2 2 2; 3 3 3]
R[:, :, 3] = [1 1 1; 2 2 2; 3 3 3]
S[:, :, 1] = [4 4 4; 4 4 4; 4 4 4]
S[:, :, 2] = [5 5 5; 5 5 5; 5 5 5]
S[:, :, 3] = [6 6 6; 6 6 6; 6 6 6]

@test meshgrid(XV, YV, XV) == (Q, R, S)


@test ndgrid(XV, YV) == ([4 4 4; 5 5 5; 6 6 6], [1 2 3; 1 2 3; 1 2 3])
@test ndgrid([XV, YV, XV]) ==   [[4, 5, 6], [1, 2, 3], [4, 5, 6]]
@test ndgrid(XV) == [4, 5, 6]

# isdefined
C = 1
@test @isdefined C
@test (@isdefined T) == false


## INPOLY
@test ClimateTools.leftorright(0.5,0.5, 1,0,1,1) == -1
@test ClimateTools.leftorright(1.5,.5, 1,0,1,1) == 1
@test ClimateTools.leftorright(1,0.5, 1,0,1,1) == 0

poly = Float64[0 0
               0 1
               1 1
               1 0
               0 0]'
p1 = [0.5, 0.5]
p2 = [0.5, 0.99]
p22 = [0.5, 1] # on edge
p23 = [0.5, 0] # on edge
p24 = [0, 0]   # on corner
p25 = [0, .4]   # on edge
p3 = [0.5, 1.1]

@test inpoly(p1, poly)
@test inpoly(p2, poly)
@test inpoly(p22, poly)
@test inpoly(p23, poly)
@test inpoly(p24, poly)
@test inpoly(p25, poly)
@test !inpoly(p3, poly)

# clockwise poly
poly = Float64[0 0
               1 0
               1 1
               0 1
               0 0]'

@test inpoly(p1, poly)
@test inpoly(p2, poly)
@test inpoly(p22, poly)
@test inpoly(p23, poly)
@test inpoly(p24, poly)
@test inpoly(p25, poly)
@test !inpoly(p3, poly)


# cross-over poly
poly = Float64[0 0
               1 0
               0 1
               1 1
               0 0]'
if VERSION >= v"0.5-"
    eval(:(@test_broken inpoly(p1, poly) )) # should be true
end
@test inpoly(p2, poly)
@test inpoly(p22, poly)
@test inpoly(p23, poly)
@test inpoly(p24, poly)
@test !inpoly(p25, poly) # different
@test !inpoly(p3, poly)


# with interior region
poly = Float64[0 0
               # interior
               0.1 0.1
               0.1 0.6
               0.6 0.6
               0.6 0.1
               0.1 0.1
               # exterior
               0 0
               0 1
               1 1
               1 0
               0 0]'
# inside interior poly: i.e. labeled as outside
@test !inpoly([0.3,0.3], poly)
@test !inpoly([0.3,0.5], poly)

poly = Float64[0 0
               # interior
               0.1 0.1
               0.1 0.6
               0.6 0.6
               # in-interior
               0.4 0.4
               0.4 0.2
               0.2 0.2
               0.2 0.4
               0.4 0.4
               # interior
               0.6 0.6
               0.6 0.1
               0.1 0.1
               # exterior
               0 0
               0 1
               1 1
               1 0
               0 0]'
# inside in-interior poly
@test inpoly([0.3,0.3], poly)
@test !inpoly([0.3,0.5], poly)

poly = Float64[0 0
               # interior
               0.1 0.1
               0.1 0.6
               0.6 0.6
               # in-interior
               0.4 0.4
               0.2 0.4
               0.2 0.2
               0.4 0.2
               0.4 0.4
               # interior
               0.6 0.6
               0.6 0.1
               0.1 0.1
               # exterior
               0 0
               0 1
               1 1
               1 0
               0 0]'
# inside in-interior poly
@test inpoly([0.3,0.3], poly)
@test !inpoly([0.3,0.5], poly)

poly = Float64[0 0
               # interior #1
               0.1 0.1
               0.1 0.6
               0.4 0.6
               0.4 0.6
               0.4 0.1
               0.1 0.1
               0 0
               # interior #2
               0.6 0.4
               0.6 0.6
               0.8 0.6
               0.8 0.4
               0.6 0.4
               0 0
               # exterior
               0 1
               1 1
               1 0
               0 0]'
@test !inpoly([0.2,0.4], poly)
@test !inpoly([0.3,0.15], poly)
@test inpoly([0.5,0.4], poly)
@test inpoly([0.5,0.2], poly)
@test !inpoly([0.7,0.5], poly)

# Test Interpolation/regridding
filename = joinpath(dirname(@__FILE__), "data", "sresa1b_ncar_ccsm3-example.nc")
C = load(filename, "tas")
# Get lat lon vector
lat = Float32.(C[1][Axis{:lat}][:])
lon = Float32.(C[1][Axis{:lon}][:])
latgrid = Float32.(C.latgrid)
longrid = Float32.(C.longrid)
# Shift longitude by 1
lon += Float32(1.0)
longrid += Float32(1.0)
axisdata = AxisArray(C[1].data, Axis{:lon}(lon), Axis{:lat}(lat), Axis{:time}(C[1][Axis{:time}][:]))
C2 = ClimGrid(axisdata, variable = "tas", longrid=longrid, latgrid=latgrid, msk=C.msk)
@test regrid(C, C2)[1].data[1, 1, 1] == 219.2400638156467
@test regrid(C, C2, min=0.0, max=0.0)[1].data[1, 1, 1] == 0.0
@test regrid(C, lon, lat)[1].data[1, 1, 1] == 219.2400638156467
@test regrid(C, lon, lat, min=0.0, max=0.0)[1].data[1, 1, 1] == 0.0

# Test applymask
# 1-D data
data = randn(3)
mask = [NaN; 1.;1.]
@test isnan(applymask(data, mask)[1])
@test applymask(data, mask)[2] == data[2]
@test applymask(data, mask)[3] == data[3]

# 2-D data
data = randn(3, 2)
mask = [[NaN; 1; 1] [1.; NaN;1.]]
@test isnan(applymask(data, mask)[1, 1]) && isnan(applymask(data, mask)[2, 2])
@test applymask(data, mask)[2, 1] == data[2, 1]
@test applymask(data, mask)[1, 2] == data[1, 2]
@test applymask(data, mask)[3, 1] == data[3, 1]
@test applymask(data, mask)[3, 2] == data[3, 2]

# 3-D data
data = randn(3, 2, 3)
mask = [[NaN; 1; 1] [1.; NaN;1.]]
@test isnan(applymask(data, mask)[1, 1, 1]) && isnan(applymask(data, mask)[2,2,1]) && isnan(applymask(data, mask)[1, 1, 2]) && isnan(applymask(data, mask)[2, 2, 2]) && isnan(applymask(data, mask)[1,1,3]) && isnan(applymask(data, mask)[2, 2, 3])

for i = 1:size(data, 1)
    @test applymask(data, mask)[2,1,i] == data[2,1,i]
    @test applymask(data, mask)[1,2,i] == data[1,2,i]
    @test applymask(data, mask)[3,1,i] == data[3,1,i]
    @test applymask(data, mask)[3,2,i] == data[3,2,i]
end

# 4-D data
data = randn(3,2,1,3)
mask = [[NaN; 1; 1] [1.; NaN;1.]]
@test isnan(applymask(data, mask)[1, 1, 1, 1]) && isnan(applymask(data, mask)[2,2,1,1]) && isnan(applymask(data, mask)[1, 1, 1,2]) && isnan(applymask(data, mask)[2, 2, 1,2]) && isnan(applymask(data, mask)[1, 1, 1,3]) && isnan(applymask(data, mask)[2, 2, 1,3])

for i = 1:size(data, 1)
    @test applymask(data, mask)[2,1,1,i] == data[2,1,1,i]
    @test applymask(data, mask)[1,2,1,i] == data[1,2,1,i]
    @test applymask(data, mask)[3,1,1,i] == data[3,1,1,i]
    @test applymask(data, mask)[3,2,1,i] == data[3,2,1,i]
end

# Test sumleapyear with StepRange{Date,Base.Dates.Day} type
d = Date(2003,1,1):Date(2008,12,31)
@test ClimateTools.sumleapyear(d) == 2

# Test timeresolution and pr_timefactor
filename = joinpath(dirname(@__FILE__), "data", "sresa1b_ncar_ccsm3-example.nc")
timevec = NetCDF.ncread(filename, "time")
@test ClimateTools.pr_timefactor(ClimateTools.timeresolution(timevec)) == 1.
@test ClimateTools.pr_timefactor("24h") == 86400.0
@test ClimateTools.pr_timefactor("12h") == 43200.0
@test ClimateTools.pr_timefactor("6h") == 21600.0
@test ClimateTools.pr_timefactor("3h") == 10800.0
@test ClimateTools.pr_timefactor("1h") == 3600.0
