#module SRIDs


# all epsg / esri comes from proj4
using Proj4

const SRID_LLA_WGS84  = :EPSG4326   # EPSG code for lon lat wgs84 (GPS) 
const SRID_ECEF_WGS84 = :EPSG4978   # EPSG code for ECEF wgs84 (GPS)


# make a function that can retrieve the authority and code from the type
function srid_params(srid::Symbol)
	srid_str = string(srid)
	authority = UTF8String(matchall(r"\D+", srid_str)[1])
	code = matchall(r"\d+", srid_str)[1]
	return (code, authority)
end

# calling this something different to not overload Proj4 stuff with generated functions
@generated function get_projection{T}(::Type{Val{T}})

	println("Gen: $(T)")

	# break it into authority and code
	(code, authority) = srid_params(T)
	code = parse(Int, code)

	check_dict = authority == "EPSG" ? Proj4.epsg : authority == "ESRI" ? Proj4.esri :error("Unknown SRID Authority")
	if !haskey(check_dict, code)
		error("Proj4 does not support: $(T)")
	end

	# add the projection info
	P = Proj4.Projection(check_dict[code])

end

#end

