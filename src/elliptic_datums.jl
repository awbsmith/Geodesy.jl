import Base: ==, show

#
# Build abstract types for elliptic datums
#
abstract SimpleEllipticDatum <: SimpleDatum
@inline style{T <: SimpleEllipticDatum}(::Type{T}) = :ELLIPTIC
(==){T1 <: SimpleEllipticDatum, T2 <: SimpleEllipticDatum}(datum_a::T1, datum_b::T2) = datum_ellipse(datum_a) == datum_ellipse(datum_b)

abstract SimpleDynEllipticDatum <: SimpleDynDatum
@inline style{T <: SimpleDynEllipticDatum}(::Type{T}) = :ELLIPTIC
(==){T1 <: SimpleDynEllipticDatum, T2 <: SimpleDynEllipticDatum}(datum_a::T1, datum_b::T2) = datum_ellipse(datum_a) == datum_ellipse(datum_b)


#
# Build the static elliptic datum
#
immutable EllipticDatum{ELLIPSE} <: SimpleEllipticDatum
    ellipse::ELLIPSE
end
@inline datum_ellipse(datum::EllipticDatum) = datum.ellipse

# instantiate the ellipse type if we know we should
@compat (::Type{EllipticDatum}){ELLIPSE <: TypedEllipse}(::Type{ELLIPSE}) = EllipticDatum(ELLIPSE())


# a type alias for unknown elliptic datums
typealias UnknownEllipticDatum EllipticDatum{UnknownEllipse}
@compat (::Type{UnknownEllipticDatum})() = EllipticDatum(UnknownEllipse())
show(io::Base.IO, ::Type{UnknownEllipticDatum}) = print(io, "UnknownEllipticDatum")

# make UnknownEllipticDatum equal to unknown datum
(==)(::UnknownEllipticDatum, ::UnknownDatum) = true
(==)(::UnknownDatum, ::UnknownEllipticDatum) = true


#
# Build the dynamic elliptic datum
#
immutable EllipticDynDatum{TIME, ELLIPSE}  <: SimpleDynEllipticDatum
    time::TIME
    ellipse::ELLIPSE
end
@inline datum_ellipse(X::EllipticDynDatum) = X.ellipse
@inline datum_time(X::EllipticDynDatum) = X.ref_time


#
# Known static elliptic datums
# Should theses be type aliases of EllipticDatum_tagged{...} instead of their own type?
#
abstract TypedEllipticDatum <: SimpleEllipticDatum
(==){T <: TypedEllipticDatum}(::T, ::T) = true
(==){T1 <: TypedEllipticDatum, T2 <: TypedEllipticDatum}(::T1, ::T2) = false
@inline datum_ellipse{T <: TypedEllipticDatum}(::T) = datum_ellipse(T)

immutable WGS84 <: TypedEllipticDatum; end
@inline datum_ellipse(::Type{WGS84}) = WGS84_Ellipse

immutable NAD27 <: TypedEllipticDatum; end # is this actually dynamic?
@inline datum_ellipse(::Type{NAD27}) = CLARKE66_Ellipse

immutable ED50 <: TypedEllipticDatum; ellipse::Void; end
@inline datum_ellipse(::Type{ED50}) = HAYFORD_Ellipse

immutable OSGB36 <: TypedEllipticDatum; end # is this actually dynamic?
@inline datum_ellipse(::Type{OSGB36}) = AIRY_Ellipse

# pretty up the display
for eType in subtypes(TypedEllipticDatum)
    str = replace(string(eType), "Geodesy.", "")
    eval(:(show(io::Base.IO, ::Type{$(eType)}) = print(io, $(str))))
end


#
# Known dynamic elliptic datums
# Should theses be type aliases of EllipticDynDatum_tagged{...} instead of their own type?
#
abstract TypedDynEllipticDatum <: SimpleDynEllipticDatum
(==){T <: TypedDynEllipticDatum}(::T, ::T) = true
(==){T1 <: TypedDynEllipticDatum, T2 <: TypedDynEllipticDatum}(::T1, ::T2) = false
@inline datum_ellipse{T <: TypedDynEllipticDatum}(::T) = datum_ellipse(T)

# Australia
immutable GDA94  <:  TypedDynEllipticDatum; end
@inline datum_ellipse(::Type{GDA94}) = GRS80_Ellipse
@inline datum_time(::Type{GDA94}) = DateTime(1994)

# Europia
immutable ETRS89  <:  TypedDynEllipticDatum; end
@inline datum_ellipse(::Type{ETRS89}) = GRS80_Ellipse
@inline datum_time(::Type{ETRS89}) = DateTime(1989)

# Americania
immutable NAD83  <:  TypedDynEllipticDatum; end
@inline datum_ellipse(::Type{NAD83}) = GRS80_Ellipse
@inline datum_time(::Type{NAD83}) = DateTime(1983)

