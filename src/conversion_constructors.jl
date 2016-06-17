# allow this function to be used as a constructor for geodesy point types
# @inline convert{Tdest <: AbstractPosition, Tsrc <: AbstractPosition}(::Type{Tdest}, X::Tsrc) = geotransform(Tdest, X)  N.B. conflicts with convert(Type{#T<:Any}, #T<:Any) at essentials.jl:59. :-(


function add_converstion_constructors(ptype, existing_types=setdiff(subtypes(AbstractPosition), [UnknownPosition, ptype]))

    # start a code block for conversion constructors
    qb = quote

        # add self conversions
        @inline convert{T <: $(ptype)}(::Type{T}, X::T) = X
        @inline convert(::Type{$(ptype)}, X::$(ptype)) = X
    end

    # add conversions to other types
    for otype in existing_types
        qn = quote
            @inline convert{Tdest <: $(otype)}(::Type{Tdest}, X::$(ptype)) = geotransform(Tdest, X)
            @inline convert{Tdest <: $(otype)}(::Type{Tdest}, X::$(ptype), crs) = geotransform(Tdest, get_crs(Tdest), X, crs)
        end
        append!(qb.args, qn.args)
    end

    #
    # add special cases here.
    #
    if (ptype != ENU)
        qn = quote

            # If the constructor has the form ENU(X, X2), and X2 is either a position or a position datum, then X2 is probably the datum for the output enu
            convert{Tdest <: ENU}(::Type{Tdest}, X::$(ptype), ref_pos::Union{PositionDatum, AbstractPosition}) = geotransform(Tdest, CRS(ENU_CS, ref_pos), X, get_crs(X))

            # If the constructor has the form LLA(X::ENU, X2), and X2 is either a positon or a position datum, then X2 is probably the datum for the input enu
            convert{oT <: $(ptype)}(::Type{oT}, X::ENU, ref_pos::Union{PositionDatum, AbstractPosition}) = geotransform(oT, get_crs(oT), X, CRS(ENU_CS, ref_pos))

        end
        append!(qb.args, qn.args)
    end
    return qb
end

#
# invoke the above
#
for ptype in setdiff(subtypes(AbstractPosition), [UnknownPosition])
    eval(add_converstion_constructors(ptype))
end

#
# can add a more general version of funny ENU constructors after the above are defined
#
# option 1
# @compat (::Type{T}){T <: ENU}(X::AbstractPosition, ref_pos::Union{PositionDatum, AbstractPosition}) = geotransform(T, CRS(ENU_CS, ref_pos), X, get_crs(X))
# option 2
@compat (::Type{T}){T <: ENU}(X::AbstractVector, ref_pos::Union{PositionDatum, AbstractPosition}) = T(X, CRS(ENU_CS(), ref_pos))            # avoid overload ambiguity
@compat (::Type{T}){T <: ENU, eT}(X::FixedVector{3,eT}, ref_pos::Union{PositionDatum, AbstractPosition}) = T(X, CRS(ENU_CS(), ref_pos))     # avoid overload ambiguity
@compat (::Type{T}){T <: ENU, eT}(X::NTuple{3,eT}, ref_pos::Union{PositionDatum, AbstractPosition}) = T(X, CRS(ENU_CS(), ref_pos))     # avoid overload ambiguity


@compat (::Type{T}){T <: ENU}(X, ref_pos::Union{PositionDatum, AbstractPosition}) = geotransform(T, CRS(ENU_CS, ref_pos), X, get_crs(X))

