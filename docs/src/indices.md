# Indices and Aggregations

All functions below accept YAXArray/Cube inputs.

## Daily Aggregations

```julia
dm = daymean(cube)
ds = daysum(cube)
```

## Yearly Aggregations

```julia
ymean = yearly_resample(cube; fct=mean)
ymax = annualmax(cube)
ymin = annualmin(cube)
ysum = annualsum(cube)
```

## Threshold and Count Indices

```julia
sdays = summerdays(cube)
fdays = frostdays(cube)
pr1 = prcp1(cube)
over = customthresover(cube, 20)
under = customthresunder(cube, 0)
```

## Ensemble Summaries

```julia
stats = ensemble_stats(cube; dim="time")
```
