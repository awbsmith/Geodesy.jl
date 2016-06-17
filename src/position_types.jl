import Base: getindex, convert, eltype, vec, ==, show

#
# Define some general convenience accessors
#

# get the coord system
get_coord_system{T <: AbstractPosition}(::Type{T}) = get_coord_system(get_crs(T))
get_coord_system(X::AbstractPosition) = get_coord_system(get_crs(X))

# get the dataum
get_datum{T <: AbstractPosition}(::Type{T}) = get_datum(get_crs(T))
get_datum(X::AbstractPosition) = get_datum(get_crs(X))

# does the type type have a crs specified in it?
defined_crs{T}(::Type{T}) = Val{false}

# defaults for abstract position (defaults should return types)
default_coord_system{T <: AbstractPosition}(::Type{T}) = UnknownCS
default_datum{T <: AbstractPosition}(::Type{T}) = UnknownDatum
default_crs{T <: AbstractPosition}(::Type{T}) = CRS{UnknownCS, UnknownDatum}


###########################
# Position type to create
###########################

const position_types = [  # point name, point cs, point fields, point datum abstract type
                        (:LLA,  LLA_CS,  (:lon, :lat, :alt), [AbstractDatum]),
                        (:LL,   LL_CS,   (:lon, :lat), [AbstractDatum]),
                        (:ECEF, ECEF_CS, (:x, :y, :z), [AbstractDatum]),
                        (:ENU,  ENU_CS,  (:east, :north, :up), [PositionDatum, AbstractPosition])  # allow position and PositionDatum's to be converted to CRS's
                        #(:UTM,  UTM_CS,  (:east, :north, :up), [AbstractDatum])
                       ]


###########################
# Generic Geodesy Position
###########################

immutable GeoPosition{CRS, N, T} <: AbstractPosition{CRS}
    _::NTuple{N, T}
    crs::CRS
end
@inline getindex(X::GeoPosition, i::Integer) = X._[i]


# get properties
defined_crs(::Type{GeoPosition}) = Val{false}
defined_crs{CRS}(::Type{GeoPosition{CRS}}) = Val{true}
defined_crs{CRS, N}(::Type{GeoPosition{CRS, N}}) = Val{true}
defined_crs{CRS, N, eT}(::Type{GeoPosition{CRS, N, eT}}) = Val{true}


get_crs_type(::Type{GeoPosition}) = Void
get_crs_type{CRS}(::Type{GeoPosition{CRS}}) = CRS
get_crs_type{CRS,N}(::Type{GeoPosition{CRS,N}}) = CRS
get_crs_type{CRS,N,eT}(::Type{GeoPosition{CRS,N,eT}}) = CRS

get_crs(X::GeoPosition) = X.crs
# get_crs{T <: GeoPosition}(::Type{T}) = get_crs_type(T)() # damn it
get_crs{CRS}(::Type{GeoPosition{CRS}}) = CRS()
get_crs{CRS,N}(::Type{GeoPosition{CRS,N}}) = CRS()
get_crs{CRS,N,eT}(::Type{GeoPosition{CRS,N,eT}}) = CRS()

# I dont want the change length from 1, so use ndims instead
ndims(::Type{GeoPosition}) = 3  # use 3 as a default
ndims{CRS, N, T}(::GeoPosition{CRS, N, T}) = N
ndims{CRS, N, T}(::Type{GeoPosition{CRS, N, T}}) = N

# return the element type
eltype{CRS, N, T}(::GeoPosition{CRS, N, T}) = T
eltype{CRS, N, T}(::Type{GeoPosition{CRS, N, T}}) = T

# make this the default default point type
default_type{CS}(::Type{CS}) = GeoPosition
default_type{CS}(::CS) = GeoPosition

# strip template info
strip_crs{T <: GeoPosition}(::Type{T}) = GeoPosition
strip_crs{T}(::Type{T}) = T # generic version, assumes there's no CRS to strip

# add constructors for two element types
@compat (::Type{T}){T <: GeoPosition}(X1, X2, crs) = T(promote(X1, X2), crs)

# add constructors for three element types
@compat (::Type{T}){T <: GeoPosition}(X1, X2, X3::Number) = T(promote(X1, X2, X3))
@compat (::Type{T}){T <: GeoPosition}(X1, X2, X3, crs) = T(promote(X1, X2, X3), crs)

@compat (::Type{GeoPosition}){N, eT}(X::NTuple{N,eT}) = GeoPosition(X, UnknownCRS())
@compat (::Type{GeoPosition{CRS}}){CRS, N, eT}(X::NTuple{N,eT}) = GeoPosition{CRS, N, eT}(X, CRS())
@compat (::Type{T}){T <: GeoPosition, N,eT, CRS}(X::NTuple{N,eT}, ::Type{CRS}) = T(X, CRS())

@compat (::Type{GeoPosition}){N, eT, CRS}(X::NTuple{N,eT}, crs::CRS) = GeoPosition{CRS, N, eT}(X, crs)
@compat (::Type{GeoPosition{CRS}}){N, eT, CRS}(X::NTuple{N,eT}, crs::CRS) = GeoPosition{CRS, N, eT}(X, crs)
@compat (::Type{GeoPosition{CRS,N}}){N, eT, CRS}(X::NTuple{N,eT}, crs::CRS) = GeoPosition{CRS, N, eT}(X, crs)

# Vector import
@compat (::Type{T}){T <: GeoPosition}(X::AbstractVector) = T((X...))              # not always type safe
@compat (::Type{T}){T <: GeoPosition}(X::AbstractVector, crs) = T((X...), crs)    # not always type safe
@compat (::Type{GeoPosition{CRS, N, eT}}){CRS, N, eT}(X::Vector{eT}) = GeoPosition{CRS, N, eT}((X...), CRS())
@compat (::Type{GeoPosition{CRS, N, eT}}){CRS, N, eT}(X::Vector{eT}, crs) = GeoPosition{CRS, N, eT}((X...), crs)


@compat (::Type{T}){T <: GeoPosition, N, eT}(X::FixedVector{N, eT}) = T((X...)::NTuple{N, eT})
@compat (::Type{T}){T <: GeoPosition, N, eT}(X::FixedVector{N, eT}, crs) = T((X...)::NTuple{N, eT}, crs)


# display
#Base.show(io::Base.IO, ::Type{GeoPosition}) = print(io, "GeoPosition{CRS, N, T}")
#Base.show{CRS, N, T}(io::Base.IO, ::Type{GeoPosition{CRS, N, T}}) = print(io, "GeoPosition{$(CRS), $(N), $(T)}")

# add some conversions
vec(X::GeoPosition) = [X._...]

#
# Allow every accessor for the default geodesy position
#
function build_all_accessors()
    defined = Dict()
    qb = quote; end
    for CS in get_CS_types()
        for (i, accessor) in enumerate(get_accessor_symbols(CS))
            # check if its already defined
            if !haskey(defined, accessor)
                qn = quote
                    @inline $(accessor)(X::GeoPosition) = X._[$(i)]
                end
                append!(qb.args, qn.args)
                defined[accessor] = i
            elseif defined[accessor] != i
                error("Conflicting storage location for accessor $(accessor) in type GeoPosition")
            end
        end
    end
    return qb
end
eval(build_all_accessors())




#########################
# Build an unknown point
#########################

immutable UnknownPosition <: AbstractPosition{UnknownCRS}
    _::Tuple{}    # TODO: keep the fields?
    crs::UnknownCRS
end

# allow construction with no inputs
@compat (::Type{UnknownPosition})() = UnknownPosition(tuple(), UnknownCRS())

defined_crs(::Type{UnknownPosition}) = Val{true}

# get properties
get_crs_type(::Type{UnknownPosition}) = UnknownCRS
get_crs(::Type{UnknownPosition}) = UnknownCRS()
get_crs(::UnknownPosition) = UnknownCRS()

# I dont want the change length from 1, so use ndims instead
ndims(::Type{UnknownPosition}) = 0
ndims(::UnknownPosition) = 0

eltype(::Type{UnknownPosition}) = Float64
eltype(::UnknownPosition) = Float64

vec(X::UnknownPosition) = Vector{Float64}[]

show(io::Base.IO, ::Type{UnknownPosition}) = print(io, "UnknownPosition")
show(io::Base.IO, ::UnknownPosition) = print(io, "UnknownPosition()")

#
# now we have this, we can make a type alias for a datum based on an UnknownPosition
#
# a type alias for unknown elliptic datums
typealias UnknownPositionDatum PositionDatum{UnknownPosition}
@compat (::Type{PositionDatum})(::Type{UnknownPosition}) = PositionDatum(UnknownPosition())  # force instantiation
show(io::Base.IO, ::Type{UnknownPositionDatum}) = print(io, "UnknownPositionDatum")

# allow it to promote to an UnknownDatum
@compat (::Type{CRS}){CS}(::Type{CS}, ::UnknownPosition) = CRS(CS, UnknownDatum()) # required
@compat (::Type{CRS})(cs, ::UnknownPosition) = CRS(cs, UnknownDatum())

# make UnknownEllipticDatum equal to unknown datum
(==)(::UnknownPositionDatum, ::UnknownDatum) = true
(==)(::UnknownDatum, ::UnknownPositionDatum) = true




#####################################
# Now build some types for commonly
# used coordinate systems
#####################################

# grab properties
function get_pos_type(sym)
    ptype = ()
    for i in 1:length(position_types)
        if sym == position_types[i][1]
            return position_types[i]
        end
    end
    return ptype
end

function BuildPointList()

	# build methods for each of the point types above
    for position_type in position_types

        # build the position type
        qb = BuildPositionType(position_type...)
        eval(qb)

        # add accessors
        qb = BuildAccessors(position_type...)
        eval(qb)

        # add constructors
        qb = BuildConstructors(position_type...)
        eval(qb)

        # add constructors that lift the datum into a CRS
        for datum_type in position_type[4]
            qb = BuildDatumLiftConstructors(position_type[1:3]..., datum_type)
            eval(qb)
        end

        # for when the CRS parameter is ommitted
        qb = BuildElementLiftConstructors(position_type...)
        eval(qb)

        # build conversions
        qb = BuildConversions(position_type...)
        eval(qb)
    end
end

#
# function to generate code to construct the point type
#
function BuildPositionType{CSType}(TypeName, ::Type{CSType}, FieldNames, args...)

    # build code for the fields
    homog_fields = [:($(field)::T) for field in FieldNames]    # homogenous field types

    qb = quote

        #
        # build it
        #
        immutable $(TypeName){CRS, T} <: AbstractPosition{CRS}
            $(homog_fields...)
            crs::CRS
        end
		@inline getindex(X::$(TypeName), i::Integer) = X.(i)

        # the number of non-crs fields
        ndims{T <: $(TypeName)}(::Type{T}) = $(length(FieldNames))
        ndims(::$(TypeName)) = $(length(FieldNames))

        # the defaults for this position type
        @inline default_coord_system{T <: $(TypeName)}(::Type{T}) = $(CSType)
        @inline default_datum(::Type{$(TypeName)}) = UnknownDatum
        @inline default_crs(::Type{$(TypeName)}) = CRS{$(CSType), UnknownDatum}

        # set this as the default point type for this CS
        default_type(::Type{$(CSType)}) = $(TypeName)
        default_type(::$(CSType)) = $(TypeName)

        # shortcuts to see if the CRS is defined in the template (there must be an easier way!)
        defined_crs(::Type{$(TypeName)}) = Val{false}
        defined_crs{CRS}(::Type{$(TypeName){CRS}}) = Val{true}
        defined_crs{CRS, eT}(::Type{$(TypeName){CRS, eT}}) = Val{true}

        # strip template info
        strip_crs{T <: $(TypeName)}(::Type{T}) = $(TypeName)

        # or add it
        #add_crs{T <: $(TypeName), CRS}(::Type{T}, ::Type{CRS}) = strip_crs(T){CRS}
        #add_crs{T <: $(TypeName), CRS}(::Type{T}, ::CRS) = strip_crs(T)

    end
end


#
# function to generate code to construct the point type
#
function BuildAccessors{CSType}(TypeName, ::Type{CSType}, FieldNames, args...)

    # CRS when the datumn is unknown
    default_CRS = :(CRS($(CSType)(), UnknownDatum())) # TODO: should this be more specific (e.g. an UnknownEllipticDatum or soemthing)?

    qb = quote

        # get properties
        get_crs_type(::Type{$(TypeName)}) = Void
        get_crs_type{CRS}(::Type{$(TypeName){CRS}}) = CRS
        get_crs_type{CRS, T}(::Type{$(TypeName){CRS, T}}) = CRS

        get_crs(X::$(TypeName)) = X.crs
        get_crs(::Type{$(TypeName)}) = $(default_CRS)
        # get_crs{T <: $(TypeName)}(::Type{T}) = get_crs_type(T)() #damn it
        get_crs{CRS}(::Type{$(TypeName){CRS}}) = CRS()
        get_crs{CRS, eT}(::Type{$(TypeName){CRS, eT}}) = CRS()

        # element type
        eltype{CRS, T}(::Type{$(TypeName){CRS,T}}) = T
        eltype{CRS, T}(::$(TypeName){CRS,T}) = T

    end

    #
    # Now add fieldaccessor methods
    #
    accessors = get_accessor_symbols(CSType)
    for (i, accessor) in enumerate(accessors)
        if i <= length(FieldNames)
            push!(qb.args,
                :(  $(accessor)(X::$(TypeName)) = X.$(FieldNames[i])   )
                 )
        else
            push!(qb.args,
                :(  $(accessor)(X::$(TypeName)) = 0.0   )  #? NaN
                 )
        end
    end
    return qb

end

#
# codegen some constructors
#
function BuildConstructors{CSType}(TypeName, ::Type{CSType}, FieldNames, args...)

    homog_fields = [:($(field)::T) for field in FieldNames]    # homogenous field types
    hetro_fields = [:($(field)) for field in FieldNames]       # heterogenous field types (seperate so we can make them numbers if desired)

    # some expressions for construction
    input_fields =  [:($(field)) for field in FieldNames]
    input_indexed = [:(X[$(i)]) for i in 1:length(FieldNames)]

    # CRS when the datumn is unknown
    default_CRS_type = :(CRS{$(CSType), UnknownDatum}) # TODO: should this be more specific (e.g. an UnknownEllipticDatum or soemthing)?
    default_CRS = :(CRS($(CSType)(), UnknownDatum()))

    qb = quote

        #
        # add constructors for it - CRS specified as in the input fields (best case)
        #
        @compat (::Type{$(TypeName)}){CRS,T}($(homog_fields...), crs::CRS) = $(TypeName){CRS,T}($(input_fields...), crs)
        @compat (::Type{$(TypeName)}){CRS}($(hetro_fields...), crs::CRS) = $(TypeName){CRS, promote_type($(input_fields...))}(promote($(input_fields...))..., crs)
        convert{CRS, T}(::Type{$(TypeName)}, X::NTuple{$(length(FieldNames)), T}, crs::CRS) = $(TypeName){CRS, T}($(input_indexed...), crs)

        # user convenience methods - dont use in implementations (These only work if the CRS can be instantiated with 0 args)
        @compat (::Type{$(TypeName)}){CRS}($(hetro_fields...), ::Type{CRS}) = $(TypeName)($(input_fields...), CRS())
        convert{T, CRS}(::Type{$(TypeName)}, X::NTuple{$(length(FieldNames)), T}, ::Type{CRS}) = $(TypeName)($(input_indexed...), CRS())

        #
        # add constructors for it - no CRS given
        #
        @compat (::Type{$(TypeName)})($(hetro_fields...)) = $(TypeName)($(input_fields...), $(default_CRS))
        convert{T}(::Type{$(TypeName)}, X::NTuple{$(length(FieldNames)), T}) = $(TypeName)($(input_indexed...), $(default_CRS))

        #
        # add constructors for it - CRS given in the type
        # (user convenience methods - dont use in implementations. These only work if the CRS can be instantiated with 0 args)
        #
        @compat (::Type{$(TypeName){CRS,T}}){CRS,T}($(homog_fields...)) = $(TypeName)($(input_fields...), CRS())
        @compat (::Type{$(TypeName){CRS}}){CRS}($(hetro_fields...)) = $(TypeName)($(input_fields...), CRS())
        convert{CRS, T}(::Type{$(TypeName){CRS}}, X::NTuple{$(length(FieldNames)), T}) = $(TypeName){CRS, T}($(input_indexed...), CRS())

        #
        # CRS specified in the type and in the output
        # ensure that they match
        #
        @compat (::Type{$(TypeName){CRS}}){CRS}($(hetro_fields...), crs::CRS) = $(TypeName)($(input_fields...), CRS())
        @compat (::Type{$(TypeName){CRS1}}){CRS1, CRS2}($(hetro_fields...), crs::CRS2) = $(TypeName)($(input_fields...), CRS1(crs))

        convert{CRS, T}(::Type{$(TypeName){CRS}}, X::NTuple{$(length(FieldNames)), T}, crs::CRS) = $(TypeName){CRS, T}($(input_indexed...), crs)
        convert{CRS1, CRS2, T}(::Type{$(TypeName){CRS1}}, X::NTuple{$(length(FieldNames)), T}, crs::CRS2) = $(TypeName){CRS, T}($(input_indexed...), CRS1(crs))

        #
        # we probably dont want Integer as the element type (its just convenient to type)
        #
        @compat (::Type{$(TypeName){CRS, T}}){CRS, T <: Integer}($(homog_fields...), crs::CRS) = $(TypeName){CRS, Float64}($([:(Float64($(field))) for field in FieldNames]...), crs)

    end
end

#
# codegen some constructors that lift the datum into a CRS (avoid these in implementations)
#
function BuildDatumLiftConstructors{CSType}(TypeName, ::Type{CSType}, FieldNames, DatumType)

    homog_fields = [:($(field)::T) for field in FieldNames]    # homogenous field types
    input_fields =  [:($(field)) for field in FieldNames]

    qb = quote

        # homogenous inputs
        @compat (::Type{$(TypeName){DATUM,T}}){DATUM <: $(DatumType), T}($(homog_fields...), datum::DATUM) =
                $(TypeName){CRS{$(CSType), DATUM}, T}($(input_fields...), CRS($(CSType)(), datum))
    end
end


#
# codegen some constructors for when the first template paramters is actually the element type
#
function BuildElementLiftConstructors{CSType}(TypeName, ::Type{CSType}, FieldNames, DatumType)

    hetero_fields = [:($(field)::Number) for field in FieldNames]    # homogenous field types
    input_fields =  [:($(field)) for field in FieldNames]

    qb = quote

        # homogenous inputs
        @compat (::Type{$(TypeName){T}}){T <: Number}($(hetero_fields...)) =
                $(TypeName)($(input_fields...), CRS($(CSType)(), UnknownDatum()))
    end
end

#
# codegen some conversions
#
function BuildConversions{CSType}(TypeName, ::Type{CSType}, FieldNames, DatumType)

    # some expressions for construction
    input_fields =  [:(X.($(field))) for field in FieldNames]
    input_indexed = [:(X[$(i)]) for i in 1:length(FieldNames)]

    typed_input_fields = [:(eT(X.$(field))) for field in FieldNames]
    untyped_input_fields = [:(X.$(field)) for field in FieldNames]

    qb = quote

        #
        # no actual converstion
        #
		convert(::Type{$(TypeName)}, X::$(TypeName)) = X

        #
        # only element type conversion
        #
		convert{CRS, eT, eT_in}(::Type{$(TypeName){CRS, eT}}, X::$(TypeName){CRS, eT_in}) =	$(TypeName){CRS, eT}($(typed_input_fields...), X.crs)

        #
        # import from a tuple
        #
        @compat (::Type{T}){T <: $(TypeName), eT}(X::NTuple{$(length(FieldNames)), eT}) = T($(input_indexed...))
        @compat (::Type{T}){T <: $(TypeName), eT}(X::NTuple{$(length(FieldNames)), eT}, crs) = T($(input_indexed...), crs)


        #
        # Standard Vectors
        #

        # export to vector
        vec(X::$(TypeName)) = [$(input_fields...)]

        # import from a vector
        @compat function (::Type{T}){T <: $(TypeName)}(X::AbstractVector)
            if (length(X) != $(length(FieldNames)))
                throw(DimensionMismatch("tried to construct a $(TypeName) with an $(length(X)) element vector"))
            end
            T($(input_indexed...))
        end

        # import from a vector and a CRS
        @compat function (::Type{T}){T <: $(TypeName)}(X::AbstractVector, crs)
            if (length(X) != $(length(FieldNames)))
                throw(DimensionMismatch("tried to construct a $(TypeName) with an $(length(X)) element vector"))
            end
            T($(input_indexed...), crs)
        end

        #
        # Fixed size arrays
        #

        # export to a FixedSizeArray
        Vec(X::$(TypeName)) = Vec($([:(X.$(field)) for field in FieldNames]...))

        # import from a Vec
        @compat (::Type{T}){T <: $(TypeName), eT}(X::FixedVector{$(length(FieldNames)), eT}) = T($(input_indexed...))
        @compat (::Type{T}){T <: $(TypeName), eT}(X::FixedVector{$(length(FieldNames)), eT}, crs) = T($(input_indexed...), crs)


		#
        # allow conversion between actual and unknown datums
        #
        convert{CS, eT, eT_in}(::Type{$(TypeName){CRS{CS, UnknownDatum}, eT}}, X::$(TypeName){CRS{CS, UnknownDatum}, eT_in}) =
            $(TypeName)($(typed_input_fields...), X.crs)

		# convert from an unknown datum to a known datum
		convert{CS, DATUM, eT, eT_in}(::Type{$(TypeName){CRS{CS, DATUM}, eT}}, X::$(TypeName){CRS{CS, UnknownDatum}, eT_in}) =
			$(TypeName)($(typed_input_fields...), CRS(get_coord_system(X.crs), DATUM()))
		convert{CS, DATUM, eT_in}(::Type{$(TypeName){CRS{CS, DATUM}}}, X::$(TypeName){CRS{CS, UnknownDatum}, eT_in}) =
			$(TypeName)($(untyped_input_fields...), CRS($(CSType), DATUM()))

		# convert from a known datum to an unknown datum
        convert{CS, DATUM, eT, eT_in}(::Type{$(TypeName){CRS{CS, UnknownDatum}, eT}}, X::$(TypeName){CRS{CS, DATUM}, eT_in}) =
			$(TypeName)($(typed_input_fields...), CRS(get_coord_system(X.crs), UnknownDatum()))
		convert{CS, DATUM, eT_in}(::Type{$(TypeName){CRS{CS, UnknownDatum}}}, X::$(TypeName){CRS{CS, DATUM}, eT_in}) =
			$(TypeName)($(untyped_input_fields...), CRS(get_coord_system(X.crs), UnknownDatum()))

    end

    return qb

end


#
# And build the point types
#
BuildPointList()


#
# Allow LL to be constructed with 3 values
#
alt_tolerance() = 1e-6   # allow this altitude to be this much from zeros


@compat @inline function (::Type{LL})(lon, lat, alt, CRS)
    @assert (alt <= alt_tolerance()) "Input altitude to the LL constructor was not zero"
    LL(lon, lat, CRS)
end

@compat @inline function (::Type{LL{T1, T2}}){T1 <: Number, T2}(lon, lat, alt::T1)
    @assert (alt <= alt_tolerance()) "Input altitude to the LL constructor was not zero"
    LL(lon, lat)
end
@compat @inline (::Type{T}){T <: LL, U}(X::NTuple{3, U}) = T(X[1], X[2])
@compat @inline (::Type{T}){T <: LL, U}(X::NTuple{3, U}, crs) = T(X[1], X[2], crs)


#
# Allow LLA to be constructed with 2 values
#
@compat @inline (::Type{T}){T <: LLA}(lon, lat) = T(lon, lat, zero(lon))
@compat @inline (::Type{T}){T <: LLA, U}(X::NTuple{2, U}) = T(X[1], X[2])
@compat @inline (::Type{T}){T <: LLA, U}(X::NTuple{2, U}, crs) = T(X[1], X[2], crs)


###############################################################################
# Method to assemble a point
# Unlike the constructor, this make decisions on whether or not to include
# the CRS in the assemble point
###############################################################################

# oT is the point to contruct, X is a tuple with the values, and crs is the output CRS
# there's no user specification on whether the user wants the CRS in the output type.
# TODO: write a function like has_size to see if the CRS is zero sized
assemble_position{pType}(::Type{pType}, X, crs) = assemble_position(pType, X, crs, defined_crs(pType))

# versions for when the output type has a CRS defined in it
function assemble_position{pType}(::Type{pType}, X, crs::CRS, defined_crs::Type{Val{true}})
    crs_type = get_crs_type(pType)
    crs = crs_type(crs)
    pType(X, crs)
end

assemble_position{pType}(::Type{pType}, X, crs, defined_crs::Type{Val{true}}) = pType(X)  # not sure what to do here, so be conservative

# versions for when the output has no CRS defiend in it
@generated function assemble_position{pType, CRSt}(::Type{pType}, X, crs::CRSt, defined_crs::Type{Val{false}})
    if (sizeof(CRSt) == 0)  # I could make a zero siz or something
        :($(pType)(X, crs))
    else
        :($(pType)(X)) # don't include the CRS by default if it has size
    end
end

