
# parameters used for distance calculations and coordinate system transforms

#
# Ellipsoid
#

immutable Ellipsoid
    a::Float64        # Semi-major axis
    b::Float64        # Semi-minor axis
    e²::Float64       # Eccentricity squared
    e′²::Float64      # Second eccentricity squared
end

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

# A few common Ellipsoids
const eWGS84      = Ellipsoid(a = "6378137.0", f_inv = "298.257223563")
const eGRS80      = Ellipsoid(a = "6378137.0", f_inv = "298.257222100882711243") # f_inv derived
const eHayford    = Ellipsoid(a = "6378388.0", f_inv = "297.0")
const eClarke1866 = Ellipsoid(a = "6378206.4",   b = "6356583.8")
const eAiry       = Ellipsoid(a = "6377563.396", b = "6356256.909")

#
# Datum.  Datums can be pretty much anything, so make them an abstract type and add to them as needed
#

abstract Datum


abstract Ellipse <: Datum  # Ellipse's are a common datum

# A few common elliptic datums
immutable WGS84 <: Ellipse end
ellipsoid(::Type{WGS84}) = eWGS84

immutable ETRS89 <: Ellipse end
ellipsoid(::Type{ETRS89}) = eGRS80

immutable NAD83 <: Ellipse end
ellipsoid(::Type{NAD83}) = eGRS80

immutable ED50 <: Ellipse end
ellipsoid(::Type{ED50}) = eHayford

immutable OSGB36 <: Ellipse end
ellipsoid(::Type{OSGB36}) = eAiry

immutable NAD27 <: Ellipse end
ellipsoid(::Type{NAD27}) = eClarke1866

# generic version
ellipsoid{T}(::Type{T}) = error("Type $(T) is not a known ellipsoid")

# get proj4 projections for each ellipse
@generated function lla_ellipse_proj{T <: Ellipse}(::Type{T})
	println("Gen: LLA $(T)")
	proj_str = @sprintf("+proj=longlat +a=%0.19f +b=%0.19f +no_defs", ellipsoid(T).a, ellipsoid(T).b)
	proj = Proj4.Projection(proj_str)
	return :($proj)
end
@generated function ecef_ellipse_proj{T <: Ellipse}(::Type{T})
	println("Gen: ECEF $(T)")
	proj_str = @sprintf("+proj=geocent +a=%0.19f +b=%0.19f +no_defs", ellipsoid(T).a, ellipsoid(T).b)
	proj = Proj4.Projection(proj_str)
	return :($proj)
end


# type for custom geoids
abstract GeoidFiles <: Datum




