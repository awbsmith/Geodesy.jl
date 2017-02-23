using FixedSizeArrays  # to do maths on points


##############################################
# A type to hold properties of this package
##############################################

# add stuff as needed
type GeodesyProperties
    geoid_dir::ASCIIString
end
geodesy_properties = GeodesyProperties(
                                       find_geoid_dir()
                                      )



######################################################
# Add types to facilitate dispatching
# the transformation calculation to other functions
######################################################


# define an abstract type to control what package handles what
abstract  AbstractPackageHandler

immutable UnknownHandler <: AbstractPackageHandler; end     # ????
immutable GeodesyHandler <: AbstractPackageHandler; end     # for point handled in this package


# and a function to return them
get_handler(X) = get_handler(typeof(X))
get_handler{T}(::Type{T}) = UnknownHandler


###############################
# Build some sort of heirarchy
# for the types
###############################

# abstract form for world coordinates
abstract  WorldPosition  <: GeodesyType

# abstract form for world surface coordinates
abstract  WorldSurfacePosition <: WorldPosition

# Heights relative to something... (basis for a compound coordinate reference system)
abstract  WorldHeight  <: WorldPosition

# abstract form for local coordinates (needs a refernce to transform to world)
abstract  LocalPosition  <: GeodesyType



# set up some helpers for the type methods
has_ellipse{T}(::Type{T}) = Val{false}                               # set true if there's an ellipse in the parameterization
has_refloc{T}(::Type{T})  = Val{false}                               # set true if there's an reference location in the parameterization
has_srid{T}(::Type{T})    = Val{false}                               # set true if there's an srid in the parameterization
has_geoid{T}(::Type{T})   = Val{false}                               # set true if there's a geoid in the parameterization
has_alt{T}(::Type{T})     = Val{false}                               # set true if there's an altitude coordinate of some kind

# and default parameters
default_params{T}(::Type{T}) = error("Default parameters not supplied for type $(T). Please overload Geodesy.default_params()")


########################
### World locations  ###
########################


# proj4 backs this and is lon lat ordering
"""
Point in Longitude-Latitude-Altitude (LLA) coordinates defined for the specified datum / ellipse.  The latitude coordinate is a geodetic latitude.

Use LLA_NULL(lon, lat, alt) if you don't want to encode the reference ellipse in the type
"""
immutable LLA{T <: AbstractDatum} <: WorldPosition
    lon::Float64
    lat::Float64
    alt::Float64
end

# useful shortcuts
typealias LLA_WGS84 LLA{WGS84}
typealias LLA_NULL LLA{UnknownDatum}

default_params{T <: LLA}(::Type{T}) = (UnknownDatum,)

# trait style functions
has_ellipse{T <: LLA}(::Type{T}) = Val{true}
get_handler{T <: LLA}(::Type{T}) = GeodesyHandler
has_alt{T <: LLA}(::Type{T}) = Val{true}




# Global cartesian coordinate system rotating with the Earth
"""
Cartesian cooridnates for a point on an ellipse

Warning: This is a Cartesian system centered on the ellipse's center and with axis direction specified by the ellipse.  This is not necessarily a "proper" ECEF frame, which would be centered on the Earth's center of mass with axis direction given by the International Reference Pole (IRP) and International Reference Meridian

         Example:
         the "eAiry" ellipse's center is not the Earth's center of mass, so converting from an eAiry based datum LLA{OSGB36} to ECEF{OSGB36} will not give a true ECEF position.  Use the SRID
         point type to get a true ECEF position if its required.


Use ECEF_NULL(lon, lat, alt) if you don't want to encode the reference ellipse in the type
"""
immutable ECEF{T <: AbstractDatum} <: WorldPosition
    x::Float64
    y::Float64
    z::Float64
end

# useful shortcuts
typealias ECEF_WGS84 ECEF{WGS84}
typealias ECEF_NULL ECEF{UnknownDatum}

default_params{T <: ECEF}(::Type{T}) = (UnknownDatum,)

# trait style functions
has_ellipse{T <: ECEF}(::Type{T}) = Val{true}
get_handler{T <: ECEF}(::Type{T}) = GeodesyHandler




##########################
### Surface locations  ###
##########################


### Point in Latitude-Longitude (LL) coordinates
# proj4 is lon lat ordering
"""
Point in Longitude-Latitude (LL) coordinates defined for the specified datum / ellipse.  The latitude coordinate is a geodetic latitude.

Assumes the height above the ellipsoid is 0
"""
immutable LL{T <: AbstractDatum} <: WorldSurfacePosition
    lon::Float64  # proj 4 is lon lat
    lat::Float64
    LL(x::Real, y::Real) = new(x, y)  # need to specify a constructor to stop the default constructor overwriting the FixedVectorNoTuple{2, Float64} constructors
end


# useful shortcuts
typealias LL_WGS84 LL{WGS84}
typealias LL_NULL LL{UnknownDatum}

default_params{T <: LL}(::Type{T}) = (UnknownDatum,)

# trait style functions
has_ellipse{T <: LL}(::Type{T}) = Val{true}
get_handler{T <: LL}(::Type{T}) = GeodesyHandler
has_alt{T <: LL}(::Type{T}) = Val{true}





##########################
### Surface heights  ###
##########################

# TODO something with this
immutable EllipHeight{T <: AbstractDatum} <: WorldHeight
    h::Float64
end
has_ellipse{T <: EllipHeight}(::Type{T}) = Val{true}  # trait style functions
default_params{T <: EllipHeight}(::Type{T}) = (UnknownDatum,)

# TODO: custom geoid heights
immutable GeoidHeight{T <: AbstractGeoid} <: WorldHeight
    h::Float64
end

default_params{T <: GeoidHeight}(::Type{T}) = (:(error("Always specify a geoid when using the GeoidHeight type")),)

has_geoid{T <: GeoidHeight}(::Type{T}) = Val{true}  # trait style functions
get_handler{T <: GeoidHeight}(::Type{T}) = Proj4Handler


###############################
### Local Coordinate Frames ###
###############################

# adding a layer of abstraction here to allow for ENU points with no LLA reference included in their template
"""
Unknown reference (allow for ENU points with no LLA reference included in their template)
"""
immutable UnknownRef <: WorldPosition end  # when we don't want to embed the reference frame in out Local coordinates
# show(io::IO, ::Type{UnknownRef}) = print(io, "???") # this is killing the code generation in type_methods.jl

# we're not code code gening methods for this type, so add this manually
get_datum(::Type{UnknownRef}) = UnknownDatum


### Point in East-North-Up (ENU) coordinates
# Local cartesian coordinate system
# Linearized about a reference point
"""
ENU{T}: East North Up point.  East and North lie in the reference ellipse's tangent plane at the reference point

Use ENU_NULL(e,n,u) if you don't want to encode the reference point in the type

# The template parameter should be a point in LL (ideally) or LLA, or an UnknownRef datatype
"""
immutable ENU{T} <: LocalPosition
    east::Float64
    north::Float64
    up::Float64
end

typealias ENU_NULL ENU{UnknownRef}

default_params{T <: ENU}(::Type{T}) = (UnknownRef, )

has_ellipse{T <: ENU}(::Type{T}) = Val{true}
has_refloc{T <: ENU}(::Type{T}) = Val{true}
get_handler{T <: ENU}(::Type{T}) = GeodesyHandler
has_alt{T <: ENU}(::Type{T}) = Val{true}


#ENU(x, y) = ENU(x, y, 0.0)

# TODO: wanted?
#immutable NED <: LocalPosition
#    north::Float64
#     east::Float64
#    down::Float64
#end
#NED(x, y) = NED(x, y, 0.0)

