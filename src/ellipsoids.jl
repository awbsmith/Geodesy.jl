
# parameters used for distance calculations and coordinate system transforms

#
# Ellipsoids 
#


# known make the ellipse type
abstract AbstractEllipse 

abstract KnownEllipse <: AbstractEllipse 

immutable UnknownEllipse <: AbstractEllipse end

# display ???'s for unkown ellipse for compactness
ellipsoid(::Type{UnknownEllipse}) = error("The ellipse is unknown")
show(io::IO, ::Type{UnknownEllipse}) = print(io, "???")



# customizable ellipse
immutable Ellipsoid              
    a::Float64                          # Semi-major axis
    b::Float64                          # Semi-minor axis
    e²::Float64                         # Eccentricity squared
    e′²::Float64                         # Second eccentricity squared
end


# build alias's for known ellipses
immutable CustomEllipse{T <: Ellipsoid} <: KnownEllipse end
ellipsoid{T <: Ellipsoid}(::Type{CustomEllipse{T}}) = T # grab the ellipse from the type


function Ellipsoid(; a::AbstractString="", b::AbstractString="", f_inv::AbstractString="")
    if isempty(a) || isempty(b) == isempty(f_inv)
        throw(ArgumentError("Specify parameter 'a' and either 'b' or 'f_inv'"))
    end
    if isempty(b)
        _ellipsoid_af(parse(BigFloat, a), parse(BigFloat, f_inv))
    else
        _ellipsoid_ab(parse(BigFloat, a), parse(BigFloat, b))
    end
end

function _ellipsoid_ab(a::BigFloat, b::BigFloat)
    e² = (a^2 - b^2) / a^2
    e′² = (a^2 - b^2) / b^2

    Ellipsoid(a, b, e², e′²)
end
function _ellipsoid_af(a::BigFloat, f_inv::BigFloat)
    b = a * (1 - inv(f_inv))
    _ellipsoid_ab(a, b)
end

#
# Predefined Ellipses
# 


# A few common Ellipsoids (TODO: steal values from Proj4!)
const eWGS84      = Ellipsoid(a = "6378137.0", f_inv = "298.257223563")
const eGRS80      = Ellipsoid(a = "6378137.0", f_inv = "298.257222100882711243") # f_inv derived
const eHayford    = Ellipsoid(a = "6378388.0", f_inv = "297.0")
const eClarke1866 = Ellipsoid(a = "6378206.4",   b = "6356583.8")
const eAiry       = Ellipsoid(a = "6377563.396", b = "6356256.909")



# build alias's for known ellipses
abstract DefinedEllipse <: KnownEllipse


# A few common ellipses
immutable WGS84_ELLIPSE <: DefinedEllipse end
ellipsoid(::Type{WGS84_ELLIPSE}) = eWGS84

immutable GRS80_ELLIPSE <: DefinedEllipse end
ellipsoid(::Type{GRS80_ELLIPSE}) = eGRS80

immutable HAYFORD_ELLIPSE <: DefinedEllipse end
ellipsoid(::Type{HAYFORD_ELLIPSE}) = eHayford

immutable AIRY_ELLIPSE <: DefinedEllipse end
ellipsoid(::Type{AIRY_ELLIPSE}) = eAiry

immutable CLARKE66_ELLIPSE <: DefinedEllipse end
ellipsoid(::Type{CLARKE66_ELLIPSE}) = eClarke1866



# build alias's for known datums 
immutable PsuedoDatum{T} <: KnownEllipse end

# And make some datums based on them
typealias WGS84 PsuedoDatum{WGS84_ELLIPSE}     
show(io::IO, ::Type{WGS84}) = print(io, "WGS84")

typealias NAD27 PsuedoDatum{CLARKE66_ELLIPSE}  # is this actually dynamic?
show(io::IO, ::Type{NAD27}) = print(io, "NAD27")

typealias ED50 PsuedoDatum{HAYFORD_ELLIPSE}
show(io::IO, ::Type{ED50}) = print(io, "ED50")

typealias OSGB36 PsuedoDatum{AIRY_ELLIPSE}             # Britania
show(io::IO, ::Type{OSGB36}) = print(io, "OSGB36")

# grab the ellipse from the type
ellipsoid{T}(::Type{PsuedoDatum{T}}) = ellipsoid(T) 


# The below are dynamic datums (a fixed point on the earth's surface doesn't move with continental drift)
# (aka the datum is dynamic allowing points to be statis)

immutable PsuedoDynDatum{T} <: KnownEllipse end

# Asutralia
typealias GDA94  PsuedoDynDatum{GRS80_ELLIPSE}
show(io::IO, ::Type{GDA94}) = print(io, "GRS80:1994")
ref_date(::Type{GDA94}) = DateTime(1994)

# Europia
typealias ETRS89 PsuedoDynDatum{GRS80_ELLIPSE}                  
show(io::IO, ::Type{ETRS89}) = print(io, "GRS80:1989")
ref_date(::Type{ETRS89}) = DateTime(1989)

# Americania
typealias NAD83 PsuedoDynDatum{GRS80_ELLIPSE}
show(io::IO, ::Type{NAD83}) = print(io, "NAD83:1983")
ref_date(::Type{NAD83}) = DateTime(1983)

# grab the ellipse from the type
ellipsoid{T}(::Type{PsuedoDynDatum{T}}) = ellipsoid(T) 




#
# type for custom geoids
#
abstract AbstractGeoid


immutable AusGeoid09 <: AbstractGeoid end
geoid_file(::Type{AusGeoid09}) = "ausgeoid09.pgm"



# dev tool, find stuff in the Proj4 dicts
# e.g find_match(Proj4.epsg, [r"proj=longlat", r"datum=WGS84"])
function find_match{T,U}(p4_dict::Dict{T,U}, exprs)
	
	if !isa(exprs, Vector)
		exprs = [exprs]
	end

	matched_key = Vector{T}(0)
	matched_str = Vector{U}(0)
	for key in keys(p4_dict)
		str = p4_dict[key]
		hasmatch = true	
		for expr in exprs
			hasmatch &= length(matchall(expr, str)) > 0
		end
		if hasmatch
			push!(matched_key, key)
			push!(matched_str, str)
		end
	end
	return [matched_key   matched_str]
end








