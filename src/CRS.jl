import Base: show, ==



# define an abstract for position now we have that
abstract AbstractPosition{CRS}

#
# CRS accessors
#

get_crs{T}(::T) = error("get_crs() not defined for type $(T)")
get_crs_type{T}(::Type{T}) = Void

#
# create an abstract type
#
abstract AbstractCRS
# show{T <: AbstractCRS}(io::Base.IO, ::Type{T}) = print(io,replace(string(T), "Geodesy.", ""))

#
# Create a CRS object
#
immutable CRS{CS,DATUM} <: AbstractCRS
    cs::CS
    datum::DATUM
end
@inline get_crs{T <: CRS}(X::T) = X

# construct an unknown one by default
@compat (::Type{CRS})() = CRS(UnknownCS(), UnknownDatum())

# now construct from other CRS types
# @compat (::Type{T}){T <: CRS}(crs::T) = crs  # overload hell
@compat (::Type{CRS{UnknownCS, UnknownDatum}})(crs) = CRS()
@compat (::Type{CRS{UnknownCS}})(crs) = CRS(UnknownCS(), get_datum(crs))
@compat (::Type{CRS{UnknownCS, DATUM}}){DATUM}(crs) = CRS{UnknownCS, DATUM}(UnknownCS(), get_datum(crs))
@compat (::Type{CRS{CS, UnknownDatum}}){CS}(crs) = CRS{CS, UnknownDatum}(get_coord_system(crs), UnknownDatum())
@compat (::Type{CRS{CS, DATUM}}){CS, DATUM}(crs) = CRS{CS, DATUM}(get_coord_system(crs), get_datum(crs))

# copy constructor allowing


show(io::Base.IO, ::Type{CRS}) = print(io, "CRS")
show{CS, DATUM}(io::Base.IO, ::Type{CRS{CS, DATUM}}) = print(io, "CRS{$(CS), $(DATUM)}")
show{CS, DATUM}(io::Base.IO, crs::CRS{CS, DATUM}) = print(io, "CRS($(crs.cs), $(crs.datum))")

#
# make some constructors that ensure things are instantiated
#
@compat (::Type{CRS{CS,DATUM}}){CS,DATUM}() = CRS(CS(), DATUM())
@compat (::Type{CRS}){CS, DATUM}(::Type{CS}, ::Type{DATUM}) = CRS(CS(), DATUM())
@compat (::Type{CRS}){DATUM}(cs, ::Type{DATUM}) = CRS(cs, DATUM())
@compat (::Type{CRS}){CS}(::Type{CS}, datum) = CRS(CS(), datum)



#
# accessors
#
@inline get_coord_system(crs::CRS) = crs.cs
@inline get_coord_system{T <: CRS}(::Type{T}) = get_coord_system_type(T)()

@inline get_coord_system_type(::Type{CRS}) = UnknownCS
@inline get_coord_system_type{CS}(::Type{CRS{CS}}) = CS()
@inline get_coord_system_type{CS, DATUM}(::Type{CRS{CS, DATUM}}) = CS()

@inline get_datum(crs::CRS) = crs.datum
@inline get_datum{T <: CRS}(::Type{T}) = get_datum_type(T)()

@inline get_datum_type(::Type{CRS}) = UnknownDatum
@inline get_datum_type{CS}(::Type{CRS{CS}}) = UnknownDatum
@inline get_datum_type{CS, DATUM}(::Type{CRS{CS, DATUM}}) = DATUM


# equality checks
(==)(crs1::CRS, crs2::CRS) = (crs1.cs == crs2.cs) && (crs1.datum == crs2.datum)


#
# exception for when you try to access these for other types
#
get_coord_system{T}(::Type{T}) = error("get_coord_system() not defined for type ::Type{$(T)}") # TODO: can I get a msg into a method error?
get_coord_system{T}(::T) = error("get_coord_system() not defined for type $(T)")

get_datum{T}(::Type{T}) = error("get_datum() not defined for type Type{$(T)}")
get_datum{T}(::T) = error("get_datum() not defined for type $(T)")

#
# special case, promote Ellipses -> EllipticDatums when constructing CRS's
#
@compat @inline (::Type{CRS}){CS}(::Type{CS}, ell::AbstractEllipse) = CRS(CS, EllipticDatum(ell)) # required
@compat @inline (::Type{CRS})(cs, ell::AbstractEllipse) = CRS(cs, EllipticDatum(ell))

@compat @inline (::Type{CRS}){CS}(::Type{CS}, ell::UnknownEllipse) = CRS(CS, UnknownDatum()) # required
@compat @inline (::Type{CRS})(cs, ell::UnknownEllipse) = CRS(cs, UnknownDatum())

#
# special case, promote Points -> PositionDatum when constructing CRS's
#
@compat @inline (::Type{CRS}){CS}(::Type{CS}, pos::AbstractPosition) = CRS(CS, PositionDatum(pos))  # required
@compat @inline (::Type{CRS})(cs, pos::AbstractPosition) = CRS(cs, PositionDatum(pos))


#
# make some type alias's for commonly used CRS's
#

typealias UnknownCRS CRS{UnknownCS, UnknownDatum}
show(io::Base.IO, ::Type{UnknownCRS}) = print(io, "UnknownCRS")
show(io::Base.IO, ::UnknownCRS) = print(io, "UnknownCRS()")

typealias LLA_CRS{DATUM} CRS{LLA_CS, DATUM}
@compat (::Type{LLA_CRS})() = CRS(LLA_CS(), UnknownDatum())
@compat (::Type{LLA_CRS})(datum) = CRS(LLA_CS(), datum)

typealias LL_CRS{DATUM} CRS{LL_CS, DATUM}
@compat (::Type{LL_CRS})() = CRS(LL_CS(), UnknownDatum())
@compat (::Type{LL_CRS})(datum) = CRS(LL_CS(), datum)

typealias ECEF_CRS{DATUM} CRS{ECEF_CS, DATUM}
@compat (::Type{ECEF_CRS})() = CRS(ECEF_CS(), UnknownDatum())
@compat (::Type{ECEF_CRS})(datum) = CRS(ECEF_CS(), datum)

typealias ENU_CRS{DATUM} CRS{ENU_CS, DATUM}
@compat (::Type{ENU_CRS})() = CRS(ENU_CS(), UnknownDatum())
@compat (::Type{ENU_CRS})(datum) = CRS(ENU_CS(), datum)

