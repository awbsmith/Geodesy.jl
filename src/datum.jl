
# parameters used for distance calculations and coordinate system transforms


#
# Datum.  It seems Datums can be pretty much anything, so make them an abstract type and subtype as needed
#

abstract Datum


# include SRID as a datum type (well it explains how to interpret a point anyway...)
include("SRIDs.jl")


#
# Time independent (static) datums.  If the datum is static then points fixed to the earths surface will change coordinates over time  (static datum -> dynamic fixed surface points)
#
abstract StaticDatum <: Datum  


#
# Time dependent (dynamic) datums.  If the datum is dynamic then points fixed to the earths surface will NOT change coordinates over time (dynamic datum -> static fixed surface points)
#
abstract DynamicDatum{T <: StaticDatum, Y} <: Datum  # Y is an Int describing the year
show{T <: DynamicDatum}(io::IO, ::Type{T}) = print(io, "$(get_params(T)[1]):$(get_params(T)[2])")



#
# Static ellipsoids Ellipsoid
#
abstract AbstractEllipse 

# known make the ellipse type
immutable EllipseDatum{T} <: StaticDatum end   # T should be an AbstractEllipse, but I can't get to the AbstractEllipse type Val{T::Ellipsoid} 


# customizable
immutable Ellipsoid                     # Val{T::Ellipsoid} <: AbstractEllipse would be really nice             
    a::Float64                          # Semi-major axis
    b::Float64                          # Semi-minor axis
    e²::Float64                         # Eccentricity squared
    e′²::Float64                         # Second eccentricity squared
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


# A few common Ellipsoids (TODO: steal values from Proj4!)
const eWGS84      = Ellipsoid(a = "6378137.0", f_inv = "298.257223563")
const eGRS80      = Ellipsoid(a = "6378137.0", f_inv = "298.257222100882711243") # f_inv derived
const eHayford    = Ellipsoid(a = "6378388.0", f_inv = "297.0")
const eClarke1866 = Ellipsoid(a = "6378206.4",   b = "6356583.8")
const eAiry       = Ellipsoid(a = "6377563.396", b = "6356256.909")


# build alias's for known ellipses
abstract KnownEllipse <: AbstractEllipse

# A few common elliptic datums
# N.B. using manual types here rather than Val{eWGS84} because its handy to have multiple names for the same ellipse
immutable WGS84_ <: KnownEllipse end
typealias WGS84 EllipseDatum{WGS84_}     # maybe WGS84 EllipseDatum{Val{eWGS84}}  would be cleaner
ellipsoid(::Type{WGS84}) = eWGS84
show(io::IO, ::Type{WGS84}) = print(io, "WGS84")

immutable GRS80_ <: KnownEllipse end
typealias GRS80 EllipseDatum{GRS80_}
ellipsoid(::Type{GRS80}) = eGRS80
show(io::IO, ::Type{GRS80}) = print(io, "GRS80")

immutable ETRS89_ <: KnownEllipse end
typealias ETRS89 EllipseDatum{ETRS89_}
ellipsoid(::Type{ETRS89}) = eGRS80
show(io::IO, ::Type{ETRS89}) = print(io, "ETRS89")

immutable NAD83_ <: KnownEllipse end
typealias NAD83 EllipseDatum{NAD83_}
ellipsoid(::Type{NAD83}) = eGRS80
show(io::IO, ::Type{NAD83}) = print(io, "NAD83")

immutable ED50_ <: KnownEllipse end
typealias ED50 EllipseDatum{ED50_}
ellipsoid(::Type{ED50}) = eHayford
show(io::IO, ::Type{ED50}) = print(io, "ED50")

immutable OSGB36_ <: KnownEllipse end
typealias OSGB36 EllipseDatum{OSGB36_}
ellipsoid(::Type{OSGB36}) = eAiry
show(io::IO, ::Type{OSGB36}) = print(io, "OSGB36")

immutable NAD27_ <: KnownEllipse end
typealias NAD27 EllipseDatum{NAD27_}
ellipsoid(::Type{NAD27}) = eClarke1866
show(io::IO, ::Type{NAD27}) = print(io, "NAD27")

# generic version of grabbing an elipse from a static datum
ellipsoid{T <: Datum}(::Type{T}) = error("Type $(T) is not a known ellipsoid")


# but allow grabbing from ellipse value types
ellipsoid{T <: Ellipsoid}(::Type{EllipseDatum{Val{T}}}) = T




# get proj4 projections for each ellipse
# using generated so they compile when T is known
@generated function lla_ellipse_proj{T}(::Type{EllipseDatum{T}})
	println("Gen: LLA $(T)")
	proj_str = @sprintf("+proj=longlat +a=%0.19f +b=%0.19f +no_defs", ellipsoid(T).a, ellipsoid(T).b)
	proj = Proj4.Projection(proj_str)
	return :($proj)
end
@generated function ecef_ellipse_proj{T}(::Type{EllipseDatum{T}})
	println("Gen: ECEF $(T)")
	proj_str = @sprintf("+proj=geocent +a=%0.19f +b=%0.19f +no_defs", ellipsoid(T).a, ellipsoid(T).b)
	proj = Proj4.Projection(proj_str)
	return :($proj)
end

#
# type for custom geoids
#


# TODO: something with this
abstract GeoidicDatum <: StaticDatum





#
# Dynamic ellipsoids as a datum
#

immutable DynamicEllipse{T <: EllipseDatum, Y} <: DynamicDatum{T, Y}  end            # Y is a year as an Int


get_params{T, Y}(::Type{DynamicEllipse{T,Y}}) = (T, Y)                 # I'd love this to be defined for the abstract DynamicDatum instead
typealias GDA94 DynamicEllipse{GRS80, 1994}  # not constistent with the way ellipse's were defined (each ellipse .  Better or worse?











