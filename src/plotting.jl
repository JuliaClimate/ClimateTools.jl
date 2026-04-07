function _plotting_extension_error(fname::Symbol)
	throw(ArgumentError(string(
		fname,
		" requires GeoMakie.jl and a Makie backend. Install GeoMakie and a backend such as CairoMakie or GLMakie, then load them before calling ",
		fname,
		"."
	)))
end

geomap(args...; kwargs...) = _plotting_extension_error(:geomap)
geomapfacet(args...; kwargs...) = _plotting_extension_error(:geomapfacet)
timeseriesplot(args...; kwargs...) = _plotting_extension_error(:timeseriesplot)
statsplot(args...; kwargs...) = _plotting_extension_error(:statsplot)

