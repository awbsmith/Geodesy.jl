import Base: ==, show

# parameters used for distance calculations and coordinate system transforms
abstract AbstractEllipse

# get the ellipse
datum_ellipse{T}(::T) = error("Can't retrieve ellipse parameters from an object of type T$(T)")
@inline datum_ellipse(X::AbstractEllipse) = X

# compare them
(==){T1 <: AbstractEllipse, T2 <: AbstractEllipse}(ell1::T1, ell2::T2) = get_parameters(ell1) == get_parameters(ell2)

#
# Ellipsoid
#
immutable Ellipsoid <: AbstractEllipse
    a::Float64                          # Semi-major axis
    b::Float64                          # Semi-minor axis
    e²::Float64                         # Eccentricity squared
    e′²::Float64                         # Second eccentricity squared
end
(==)(e1::Ellipsoid, e2::Ellipsoid) = (e1.a == e2.a) && (e1.b == e2.b)

# add constructors
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

# A few common Ellipsoids (TODO: check against Proj4!)
const eWGS84      = Ellipsoid(a = "6378137.0", f_inv = "298.257223563")
const eGRS80      = Ellipsoid(a = "6378137.0", f_inv = "298.257222100882711243") # f_inv derived
const eHayford    = Ellipsoid(a = "6378388.0", f_inv = "297.0")
const eClarke1866 = Ellipsoid(a = "6378206.4",   b = "6356583.8")
const eAiry       = Ellipsoid(a = "6377563.396", b = "6356256.909")

#
# Define tagged ellipses
#
abstract TypedEllipse <: AbstractEllipse
@inline get_parameters{T <: TypedEllipse}(::T) = get_parameters(T)

# overload the previously defined equality
@inline (==){T <: TypedEllipse}(ell_a::T, ell_b::T) = true
@inline (==){T1 <: TypedEllipse, T2 <: TypedEllipse}(ell_a::T1, ell_b::T2) = false

immutable UnknownEllipse <: TypedEllipse; end
get_parameters(::Type{UnknownEllipse}) = error("Parameters can't be retrieved from type UnknownEllipse")

immutable WGS84_Ellipse <: TypedEllipse end
get_parameters(::Type{WGS84_Ellipse}) = eWGS84

immutable GRS80_Ellipse <: TypedEllipse end
get_parameters(::Type{GRS80_Ellipse}) = eGRS80

immutable HAYFORD_Ellipse <: TypedEllipse end
get_parameters(::Type{HAYFORD_Ellipse}) = eHayford

immutable AIRY_Ellipse <: TypedEllipse end
get_parameters(::Type{AIRY_Ellipse}) = eAiry

immutable CLARKE66_Ellipse <: TypedEllipse end
get_parameters(::Type{CLARKE66_Ellipse}) = eClarke1866

