#
# Datums 
# note that geodesy itself doesn't understand datums, it just stops you doing silly things with them

abstract AbstractDatum

# start a get_datum function to retrieve the datum where possible
get_datum{T <: AbstractDatum}(::Type{T}) = T

#
# Unknown
#
immutable UnknownDatum <: AbstractDatum end
get_datum(::Type{AbstractDatum}) = UnknownDatum  # make the abstract datum return this


# display ???'s for unkown ellipse for compactness
ellipsoid(::Type{UnknownDatum}) = error("The ellipse is unknown")
# show(io::IO, ::Type{UnknownDatum}) = print(io, "???") # this is killing the code generation in type_methods.jl


#
# Known
#
abstract KnownDatum <: AbstractDatum




# build alias's for known datums 
"""
When using a "static" datum fixed points on the Earth's surface move over time because of effects like continental drift (the datum is static but the Earth is not)
"""
abstract StaticDatum <: KnownDatum 

# And make some datums based on them
immutable WGS84 <: StaticDatum; end
ellipsoid(::Type{WGS84}) = ellipsoid(WGS84_ELLIPSE)
show(io::IO, ::Type{WGS84}) = print(io, "WGS84")

immutable NAD27 <: StaticDatum; end # is this actually dynamic?
ellipsoid(NAD27) = ellipsoid(CLARKE66_ELLIPSE)
show(io::IO, ::Type{NAD27}) = print(io, "NAD27")

immutable ED50 <: StaticDatum; end 
ellipsoid(::Type{ED50}) = ellipsoid(HAYFORD_ELLIPSE)
show(io::IO, ::Type{ED50}) = print(io, "ED50")

immutable OSGB36 <: StaticDatum; end # is this actually dynamic?
ellipsoid(::Type{OSGB36}) = ellipsoid(AIRY_ELLIPSE)
show(io::IO, ::Type{OSGB36}) = print(io, "OSGB36")


# The below are dynamic datums (a fixed point on the earth's surface doesn't move with continental drift)
# (aka the datum is dynamic allowing points to be static)
"""
When using a "dynamic" datum fixed points on the Earth's surface have constant coordinates (static) over time, because the datum changes (dynamic) over time to account for effects like continental drift.  
Any accurate transformation between static and dynamic datums should include a time
"""
abstract DynDatum <: KnownDatum 

# Australia
immutable GDA94  <:  DynDatum; end
ellipsoid(::Type{GDA94}) = ellipsoid(GRS80_ELLIPSE)
ref_date(::Type{GDA94}) = DateTime(1994)
show(io::IO, ::Type{GDA94}) = print(io, "GDA94")

# Europia
immutable ETRS89  <:  DynDatum; end
ellipsoid(::Type{ETRS89}) = ellipsoid(GRS80_ELLIPSE)
ref_date(::Type{ETRS89}) = DateTime(1989)
show(io::IO, ::Type{ETRS89}) = print(io, "ETRS89")


# Americania
immutable NAD83  <:  DynDatum; end
ellipsoid(::Type{NAD83}) = ellipsoid(GRS80_ELLIPSE)
show(io::IO, ::Type{NAD83}) = print(io, "NAD83")
ref_date(::Type{NAD83}) = DateTime(1983)


#
# type for custom geoids
#
abstract AbstractGeoid

"""
Function to set / get the geoid directory
"""
function set_geoid_dir(path::AbstractString)
    Geodesy.geodesy_properties.geoid_dir = path
end
get_geoid_dir() = Geodesy.geodesy_properties.geoid_dir

# dont know the Geoid
immutable UnknownGeoid <: AbstractGeoid end
geoid_file(::Type{UnknownGeoid}) = error("Unsure which Geoid to use for an unknown Geoid")

# Known Geoids
abstract KnownGeoid <: AbstractGeoid

# Australia
immutable AusGeoid09 <: KnownGeoid end
geoid_file(::Type{AusGeoid09}) = "ausgeoid09.gtx"





# function to get a list of all datums
function get_datums(datum::DataType=AbstractDatum, super::ASCIIString="", out=Vector{ASCIIString}(0))
    sub_datums = subtypes(datum)
    if (length(sub_datums) > 0)
        for sd in sub_datums
            out = get_datums(sd, datum==AbstractDatum ? "" : " (" * string(datum) * ")", out)
        end
    else
        push!(out, replace(string(datum) *  super  , "Geodesy.", ""))
    end
    return out
end







#####################################################################
# dev tool, find stuff in the Proj4 dicts
# e.g find_match(Proj4.epsg, [r"proj=longlat", r"datum=WGS84"])
#####################################################################
function find_match{T,U}(p4_dict::Dict{T,U}, exprs)
    
    if !isa(exprs, Vector)
        exprs = [exprs]
    end

    matched_key = Vector{T}(0)
    matched_str = Vector{U}(0)
    for key in keys(p4_dict)
        str = p4_dict[key]
        hasmatch = true    
        for expr in exprs
            hasmatch &= length(matchall(expr, str)) > 0
        end
        if hasmatch
            push!(matched_key, key)
            push!(matched_str, str)
        end
    end
    return [matched_key   matched_str]
end
