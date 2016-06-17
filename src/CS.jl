import Base: convert, show

# define different coordinate systems
abstract AbstractCS

# when we don't know
abstract AbstractUnknownCS <: AbstractCS

#
# function to get the default point type for this CS
#
default_type{T <: AbstractCS}(::Type{T}) = error("No default point type for coordinate system: $(T)")

# define some traist for the coordinate system
cs_prop_fcns = [:is_cart, :is_polar, :is_enu, :is_geodetic, :is_timed]

# template cs properaties
function create_cs_templates()
    qb = quote; end
    for fcn in cs_prop_fcns
        push!(qb.args, :(@inline $(fcn){T}(::Type{T}) = Val{false}))
        push!(qb.args, :(@inline $(fcn){T}(::T) = $(fcn)(T)))
    end
    return qb
end
eval(create_cs_templates())

#
# Unknown coordinate systems
#
immutable UnknownCS <: AbstractUnknownCS; end

# always allow conversion to an UnknownCS
convert(::Type{UnknownCS}, cs::AbstractCS) = UnknownCS() # make this anything, or just AbstractCS?


#
# Cartesian variants
#
abstract AbstractCartCS <: AbstractCS
@inline is_cart{T <: AbstractCartCS}(::Type{T}) = Val{true}

# concrete versions
immutable ECEF_CS <: AbstractCartCS; end
@inline is_geodetic(::Type{ECEF_CS}) = Val{true}

immutable ENU_CS <: AbstractCartCS; end
@inline is_geodetic(::Type{ENU_CS}) = Val{true}


#
# Polar styles
#
abstract AbstractPolarCS <: AbstractCS
@inline is_polar{T <: AbstractPolarCS}(::Type{T}) = Val{true}

# concrete versions
immutable LL_CS <: AbstractPolarCS; end
@inline is_geodetic(::Type{LL_CS}) = Val{true}

immutable LLA_CS <: AbstractPolarCS; end
@inline is_geodetic(::Type{LLA_CS}) = Val{true}

#
# Others
#
abstract AbstractFauCartCS <: AbstractCS      # this is a lower requirement the Cartesian (i.e. not quite Cartesian)

immutable UTM_CS <: AbstractFauCartCS; end
@inline is_geodetic(::Type{UTM_CS}) = Val{true}


#
# now make a function to create the expected accessors for each CS type
#
function get_accessor_symbols{csType}(::Type{csType})
    if (is_cart(csType) == Val{true}) && (is_geodetic(csType) == Val{true})
        list = [:get_x, :get_y, :get_z]
    elseif (is_polar(csType) == Val{true}) && (is_geodetic(csType) == Val{true})
        list = [:get_lon, :get_lat, :get_alt]
    elseif (is_enu(csType) == Val{true})
        list = [:get_east, :get_north, :get_up]
    else
        list = []
    end
    if (is_timed(csType) == Val{true})
        push!(list, :get_time)
    end
    return list
end


#
# assess equality of coordinmate systems
#

@inline function combine_cs{CS}(cs1::CS, cs2::CS)
    (cs1 == cs2) || error("The coordinates systems $(cs1) and $(cs2) don't match")
    return cs1
end

@inline combine_cs(cs1::AbstractUnknownCS, ::AbstractUnknownCS) = UnknownCS()
@inline combine_cs(cs1::AbstractUnknownCS, cs2) = cs2
@inline combine_cs(cs1, cs2::AbstractUnknownCS) = cs1
function combine_cs{CS}(cs1::Type{CS}, cs2::Type{CS})
    cs1 == cs2 || error("coordinate system $(cs1) != coordinate system $(cs1)")
    return cs1
end
@inline combine_cs{CS1, CS2}(cs1::CS1, cs2::CS2) = error("coordinate system $(cs1) != coordinate system $(cs1)")

# function to return a list of all concrete CS types
function get_CS_types()
    bucket = subtypes(AbstractCS)
    list = Vector{Any}(0)
    while (length(bucket) > 0)
        cs = pop!(bucket)
        if isleaftype(cs)
            push!(list, cs)
        end
        append!(bucket, subtypes(cs))
    end
    return list
end

