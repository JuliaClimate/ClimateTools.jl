const _GEOMAKIE_PKGID = Base.PkgId(Base.UUID("db073c08-6b98-4ee5-b6a4-5efafb3259c6"), "GeoMakie")

function _ensure_plotting_extension_loaded()
	ext = Base.get_extension(@__MODULE__, :ClimateToolsGeoMakieExt)
	ext !== nothing && return true

	try
		Base.require(_GEOMAKIE_PKGID)
	catch
		return false
	end

	return Base.get_extension(@__MODULE__, :ClimateToolsGeoMakieExt) !== nothing
end

function _plotting_extension_error(fname::Symbol)
	throw(ArgumentError(string(
		fname,
		" requires GeoMakie.jl and a Makie backend. Install GeoMakie and a backend such as CairoMakie or GLMakie, then load them before calling ",
		fname,
		"."
	)))
end

function geomap(args...; kwargs...)
	_ensure_plotting_extension_loaded() && return Base.invokelatest(geomap, args...; kwargs...)
	_plotting_extension_error(:geomap)
end

function geomapfacet(args...; kwargs...)
	_ensure_plotting_extension_loaded() && return Base.invokelatest(geomapfacet, args...; kwargs...)
	_plotting_extension_error(:geomapfacet)
end

function timeseriesplot(args...; kwargs...)
	_ensure_plotting_extension_loaded() && return Base.invokelatest(timeseriesplot, args...; kwargs...)
	_plotting_extension_error(:timeseriesplot)
end

function statsplot(args...; kwargs...)
	_ensure_plotting_extension_loaded() && return Base.invokelatest(statsplot, args...; kwargs...)
	_plotting_extension_error(:statsplot)
end

function robustnessmap(args...; kwargs...)
	_ensure_plotting_extension_loaded() && return Base.invokelatest(robustnessmap, args...; kwargs...)
	_plotting_extension_error(:robustnessmap)
end

