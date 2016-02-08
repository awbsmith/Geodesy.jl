using FixedSizeArrays  # to do maths on points

using Proj4

###############################
# Build some sort of heirarchy
# for the types
###############################

# abstract form for world coordinates
abstract  WorldPosition <: FixedVectorNoTuple{3, Float64}

# abstract form for world surface coordinates
abstract  WorldSurfacePosition <: FixedVectorNoTuple{2, Float64}

# abstract form for local coordinates (needs a refernce to transform to world)
abstract  LocalPosition <: FixedVectorNoTuple{3, Float64}

# Heights relative to something...
abstract  WorldHeight

# include SRID info
include("SRIDs.jl")

##########################
### Surface locations  ###
##########################


### Point in Latitude-Longitude (LL) coordinates
# proj4 is lon lat ordering
immutable LL{T <: Datum} <: WorldSurfacePosition
	lon::Float64  # proj 4 is lon lat    
	lat::Float64
	LL(x::Real, y::Real) = new(x, y)  # need to specify a constructor to stop the default constructor overwriting the FixedVectorNoTuple{2, Float64} constructors
end
#call{T}(::Type{LL{T}}, a::Tuple, b...) = call(FixedVectorNoTuple, a, b...)  # restore the FSA style constuctor

# get the reference ellipsoid
ellipsoid{T}(::Type{LL{T}}) = ellipsoid(T)


##########################
### Surface heights  ###
##########################

immutable EllipHeight{T <: Datum} <: WorldHeight
	alt::Float64
end

# custom geoids requiring a "pgm" file
immutable GeoidHeight{T <: GeoidFiles} <: WorldHeight
	alt::Float64
end


########################
### World locations  ###
########################

# proj4 is lon lat ordering
### Point in Latitude-Longitude-Altitude (LLA) coordinates
immutable LLA{T <: Datum} <: WorldPosition
	lon::Float64    
	lat::Float64
    alt::Float64
end
Base.call{T}(::Type{LLA{T}}, lon::Real, lat::Real) = LLA{T}(lon, lat, NaN)
Base.call{T}(::Type{LLA{T}}; lat::Real=NaN, lon::Real=NaN, h::Real=h::Real) = LLA{T}(lon, lat, h)
ellipsoid{T}(::Type{LLA{T}}) = ellipsoid(T)  # reference ellipsoid for the position

# typealias common usage cases here
typealias LLA_WGS84 LLA{WGS84_ELLIPSE}

### Point in Earth-Centered-Earth-Fixed (ECEF) coordinates
# Global cartesian coordinate system rotating with the Earth
immutable ECEF{T <: Datum} <: WorldPosition 
    x::Float64
    y::Float64
    z::Float64
end
typealias ECEF_WGS84 ECEF{WGS84_ELLIPSE}


### SRID based points (converions will use Proj4)
### N.B, for lat lon style SRIDs,  x -> lon, y -> lat (or getlat() and getlon())
immutable SRID{T} <: WorldPosition
   	x::Float64
    y::Float64
	z::Float64
end
call{T}(::Type{SRID{T}}, x::Real, y::Real) = call(SRID{T}, x, y, NaN)  # Using NaN to check in case height is important for a transformation (it will pollute)

# convenience for getting the srid from an SRID point
srid_string{T}(X::SRID{T}) = string(T)
srid_params{T}(X::SRID{T}) = srid_params(T)










###############################
### Local Coordinate Frames ###
###############################


### Point in East-North-Up (ENU) coordinates
# Local cartesian coordinate system
# Linearized about a reference point
immutable ENU     <: LocalPosition
    east::Float64
    north::Float64
    up::Float64
end
#ENU(x, y) = ENU(x, y, 0.0)


immutable NED <: LocalPosition
    north::Float64
	east::Float64
    down::Float64
end
#NED(x, y) = NED(x, y, 0.0)


### XYZ
# Helper for creating other point types in generic code
# e.g. myfunc{T <: Union(ENU, LLA)}(...) = (x, y = ...; T(XY(x, y)))
#=
type XYZ
    x::Float64
    y::Float64
    z::Float64
end
XY(x, y) = XYZ(x, y, 0.0)

Base.call{T}(::Type{LL{T}}, xyz::XYZ) = LL{T}(xyz.y, xyz.x)
Base.call{T}(::Type{LLA{T}}, xyz::XYZ) = LLA{T}(xyz.y, xyz.x, xyz.z)
ENU(xyz::XYZ) = ENU(xyz.x, xyz.y, xyz.z)
=#

### get*
# Point translators
getX(ll::LL) = ll.lon
getY(ll::LL) = ll.lat

getX(lla::LLA) = lla.lon
getY(lla::LLA) = lla.lat
getZ(lla::LLA) = lla.alt

getX(enu::ENU) = enu.east
getY(enu::ENU) = enu.north
getZ(enu::ENU) = enu.up

getlat(X::SRID) = X.y
getlon(X::SRID) = X.x




