import Base: show, convert

abstract AbstractSRID <: AbstractCRS

# a data version
immutable SRID  <:  AbstractSRID
    auth::Symbol
    code::Int
end

# get the template components
get_auth(srid::SRID) = srid.auth
get_code(srid::SRID) = srid.code
@inline Base.show(io::Base.IO, ::SRID) = print(io, "SRID(:$(srid.auth), $(srid.code))")

# an unkown version
immutable UnknownSRID  <:  AbstractSRID; end
get_auth(srid::UnknownSRID) = error("The SRID's authority is unknown")
get_code(srid::UnknownSRID) = error("The SRID's code is unknown")
get_crs(srid::UnknownSRID) = UnknownCRS()

# also use a pure template version
immutable SRIDt{AUTH, CODE}  <:  AbstractSRID
    auth::Void # TODO: keep these null fields?
    code::Void
end
@compat (::Type{SRIDt{AUTH, CODE}}){AUTH, CODE}() = SRIDt{AUTH, CODE}(nothing, nothing)
@compat (::Type{SRIDt})(auth, code) = SRIDt{auth, code}(nothing, nothing)

# get the template components
get_auth{AUTH, CODE}(srid::SRIDt{AUTH, CODE}) = AUTH
get_code{AUTH, CODE}(srid::SRIDt{AUTH, CODE}) = CODE
# @inline Base.show{AUTH, CODE}(io::Base.IO, ::SRIDt{AUTH, CODE}) = print(io, "SRIDt{:$(AUTH), $(CODE)}")

# conversions
convert(::Type{SRIDt}, srid::SRID) = SRIDt{get_auth(srid), get_code(srid)}(nothing, nothing)
convert(::Type{SRID}, srid::SRIDt) = SRID(get_auth(srid), get_code(srid))


# method to relate it to a CRS
@inline get_coord_system(srid::AbstractSRID) = get_coord_system(get_crs(srid)) # N.B. this will be a DataType
@inline get_datum(srid::AbstractSRID) = get_datum(get_crs(srid))               # N.B. this will be a DataType

