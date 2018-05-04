# Maps

## Mapping the ClimGrid type

Mapping climate information can be done by using `mapclimgrid`:
```julia
mapclimgrid(C::ClimGrid; region = "World")
```

Which should return

<p align="center">
  <img src="https://cloud.githubusercontent.com/assets/3630311/23712122/e97bd322-03ef-11e7-93da-749c961c4070.png?raw=true" width="771" height="388" alt="Precipitation example"/>
</p>

Note that if the `ClimGrid` data structure has 3 dimensions (time x longitude x latitude) the `mapclimgrid` function makes a time-average (i.e. climatological mean). Right now, options are available for region: `World`, `Canada`, `Quebec` and the default `auto` which use the maximum and minimum of the lat-long coordinates inside the `ClimGrid` structure. The user can also provide a polygon(s) and the `mapclimgrid` function will clip the grid points outside the specified polygon. Another option is to provide a mask (with dimensions identical to the spatial dimension of the `ClimGrid` data) which contains `NaN` and `1.0` and the data inside the `ClimGrid` struct will be clipped with the mask. Other regions will be added in the future, as well as the option to send a custom region defined by a lat-lon box.

In a future release, the user will have the option to specify his own time period (e.g. plotting the time-average of a given month and year, as opposed to the time-average of the whole `ClimGrid` structure).
