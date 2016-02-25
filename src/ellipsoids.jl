
# parameters used for distance calculations and coordinate system transforms

#
# Ellipsoids 
#
immutable Ellipsoid              
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


#
# Predefined Ellipses
# 

# A few common Ellipsoids (TODO: steal values from Proj4!)
const eWGS84      = Ellipsoid(a = "6378137.0", f_inv = "298.257223563")
const eGRS80      = Ellipsoid(a = "6378137.0", f_inv = "298.257222100882711243") # f_inv derived
const eHayford    = Ellipsoid(a = "6378388.0", f_inv = "297.0")
const eClarke1866 = Ellipsoid(a = "6378206.4",   b = "6356583.8")
const eAiry       = Ellipsoid(a = "6377563.396", b = "6356256.909")



#
# Generate Ellipse types
# 


# build known ellipses
abstract Ellipse <: KnownDatum

# custom type
immutable CustomEllipse{T <: Ellipsoid} <: Ellipse end
ellipsoid{T <: Ellipsoid}(::Type{CustomEllipse{T}}) = T # grab the ellipse from the type

# A few common ellipses
immutable WGS84_ELLIPSE <: Ellipse end
ellipsoid(::Type{WGS84_ELLIPSE}) = eWGS84

immutable GRS80_ELLIPSE <: Ellipse end
ellipsoid(::Type{GRS80_ELLIPSE}) = eGRS80

immutable HAYFORD_ELLIPSE <: Ellipse end
ellipsoid(::Type{HAYFORD_ELLIPSE}) = eHayford

immutable AIRY_ELLIPSE <: Ellipse end
ellipsoid(::Type{AIRY_ELLIPSE}) = eAiry

immutable CLARKE66_ELLIPSE <: Ellipse end
ellipsoid(::Type{CLARKE66_ELLIPSE}) = eClarke1866












