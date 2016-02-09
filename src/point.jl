using FixedSizeArrays  # to do maths on points

using Proj4

###############################
# Build some sort of heirarchy
# for the types
###############################

# abstract form for world coordinates
abstract  WorldPosition  <: FixedVectorNoTuple{3, Float64}  

# abstract form for world surface coordinates
abstract  WorldSurfacePosition  <: FixedVectorNoTuple{2, Float64}

# abstract form for local coordinates (needs a refernce to transform to world)
abstract  LocalPosition  <: FixedVectorNoTuple{3, Float64}

# Heights relative to something...
abstract  WorldHeight <: Real

# adding a layer of abstraction here to allow for ENU points with no LLA reference included in their template
abstract  AbstractLLA <: WorldPosition

# include SRID info
include("SRIDs.jl")


########################
### World locations  ###
########################

# proj4 is lon lat ordering
### Point in Latitude-Longitude-Altitude (LLA) coordinates
immutable LLA{T <: Datum} <: AbstractLLA
	lon::Float64    
	lat::Float64
    alt::Float64
end
Base.call{T}(::Type{LLA{T}}, lon::Real, lat::Real) = LLA{T}(lon, lat, NaN)
Base.call{T}(::Type{LLA{T}}; lat::Real=NaN, lon::Real=NaN, h::Real=h::Real) = LLA{T}(lon, lat, h)




### Point in Earth-Centered-Earth-Fixed (ECEF) coordinates
# Global cartesian coordinate system rotating with the Earth
# These are parameterized because the Earth's center of mass isn't exactly known
immutable ECEF{T <: Datum} <: WorldPosition 
    x::Float64
    y::Float64
    z::Float64
end

#  common usage typealias here
typealias ECEF_WGS84 ECEF{WGS84}
typealias LLA_WGS84 LLA{WGS84}


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

# get the reference ellipsoid
ellipsoid{T}(::Type{LL{T}}) = ellipsoid(T)


##########################
### Surface heights  ###
##########################

# TODO something with this
immutable EllipHeight{T <: Datum} <: WorldHeight
	h::Float64
end

# TODO: custom geoid heights
#immutable GeoidHeight{T <: GeoidFiles} <: WorldHeight
#	h::Float64
#end


###############################
### Local Coordinate Frames ###
###############################

type LLA_Null <: AbstractLLA end  # when we don't want to embed the reference frame in out Local coordinates

### Point in East-North-Up (ENU) coordinates
# Local cartesian coordinate system
# Linearized about a reference point
immutable ENU{T} <: LocalPosition
    east::Float64
    north::Float64
    up::Float64
end
LLA_ref{T}(::ENU{T}) = T
call(::Type{ENU}, e::Real, n::Real, u::Real) = ENU{LLA_Null}(e,n,u)  # allow default constructuction with no reference LLA
#ENU(x, y) = ENU(x, y, 0.0)


# TODO: wanted?
#immutable NED <: LocalPosition
#    north::Float64
#	 east::Float64
#    down::Float64
#end
#NED(x, y) = NED(x, y, 0.0)

#=
Base.call{T}(::Type{LL{T}}, xyz::XYZ) = LL{T}(xyz.y, xyz.x)
Base.call{T}(::Type{LLA{T}}, xyz::XYZ) = LLA{T}(xyz.y, xyz.x, xyz.z)
ENU(xyz::XYZ) = ENU(xyz.x, xyz.y, xyz.z)
=#

# retrieve datums and ellipsoids
ellipsoid{T <: Ellipse}(::Union{LLA{T}, LL{T}, ECEF{T}}) = ellipsoid(T)  # reference ellipsoid for the position
datum{T <: Datum}(::Union{LLA{T}, LL{T}, ECEF{T}}) = T                   # reference datum for the position

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

getlon(X::SRID) = X.x
getlat(X::SRID) = X.y





