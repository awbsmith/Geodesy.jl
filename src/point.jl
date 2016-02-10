using FixedSizeArrays  # to do maths on points
using Proj4

###############################
# Build some sort of heirarchy
# for the types
###############################

# abstract form for world coordinates
abstract  WorldPosition  <: FixedVectorNoTuple{3, Float64}  

# abstract form for static world coordinates, where a fixed point on the Earth's surface has the same coordinate every year (adjusts for continental drift etc).  
# E.g. ITRS syles
# Transformations from StaticWorldPosition <-> WorldPosition require a date
abstract  StaticWorldPosition  <: FixedVectorNoTuple{3, Float64}

# abstract form for world surface coordinates
abstract  WorldSurfacePosition  <: FixedVectorNoTuple{2, Float64}

# abstract form for local coordinates (needs a refernce to transform to world)
abstract  LocalPosition  <: FixedVectorNoTuple{3, Float64}

# Heights relative to something...
abstract  WorldHeight <: Real




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



### Point in Earth-Centered-Earth-Fixed (ECEF) coordinates
# Global cartesian coordinate system rotating with the Earth
# These are parameterized by a datum because the Earth's center of mass isn't exactly known
immutable ECEF{T <: Datum} <: WorldPosition
    x::Float64
    y::Float64
    z::Float64
end


### SRID based points (converions will use Proj4)
### N.B, for lat lon style SRIDs,  x -> lon, y -> lat (or getlat() and getlon())
### N.B, for utm style SRIDs,  x -> east, y -> north, z -> up (or geteast() getnorth() getup())
immutable SRID_Pos{T <: SRID} <: WorldPosition
   	x::Float64
    y::Float64
	z::Float64
end
call{T}(::Type{SRID_Pos{T}}, x::Real, y::Real) = call(SRID_Pos{T}, x, y, NaN)  # Using NaN to check in case height is important for a transformation (it will pollute)

# convenience for getting the srid from an SRID point
get_srid{T}(X::SRID_Pos{T}) = T


###############################
### Static World locations  ###
###############################

immutable StaticLLA{T <: DynamicDatum} <: StaticWorldPosition  # Y is a year as an Int
	lon::Float64    
	lat::Float64
    alt::Float64
end



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

# adding a layer of abstraction here to allow for ENU points with no LLA reference included in their template
immutable NullPos <: WorldPosition end  # when we don't want to embed the reference frame in out Local coordinates

### Point in East-North-Up (ENU) coordinates
# Local cartesian coordinate system
# Linearized about a reference point
immutable ENU{T} <: LocalPosition   # T should be either NullLLA or something along the lines of a Val{WorldPosition}
    east::Float64
    north::Float64
    up::Float64
end
call(::Type{ENU}, e::Real, n::Real, u::Real) = ENU{NullPos}(e,n,u)  								# allow default constructuction with no reference LLA
convert{T}(::Type{ENU{NullPos}}, enu::Type{ENU{T}}) = ENU{NullPos}(enu.east, enu.north, enu.up)		# allow stripping the reference position out of the template
#ENU(x, y) = ENU(x, y, 0.0)



# TODO: wanted?
#immutable NED <: LocalPosition
#    north::Float64
#	 east::Float64
#    down::Float64
#end
#NED(x, y) = NED(x, y, 0.0)

##
#  common usage typealias here
##
typealias ECEF_WGS84 ECEF{WGS84}
typealias LLA_WGS84 LLA{WGS84}
typealias LLA_GDA94 StaticLLA{GDA94}

# get srids known point / dataum combo
get_srid(::Type(LLA_WGS84)) = SRID{:EPSG, 4326}   # EPSG code for lon lat wgs84 (GPS).  This may be updated later
get_srid(::Type(ECEF_WGS84)) = SRID{:EPSG, 4978} 


#=
Base.call{T}(::Type{LL{T}}, xyz::XYZ) = LL{T}(xyz.y, xyz.x)
Base.call{T}(::Type{LLA{T}}, xyz::XYZ) = LLA{T}(xyz.y, xyz.x, xyz.z)
ENU(xyz::XYZ) = ENU(xyz.x, xyz.y, xyz.z)
=#

# retrieve datums and ellipsoids
ellipsoid{T <: Ellipse}(::Union{LLA{T}, LL{T}, ECEF{T}}) = ellipsoid(T)  # reference ellipsoid for the position
datum{T <: Datum}(::Union{LLA{T}, LL{T}, ECEF{T}, StaticLLA{T}}) = T                   # reference datum for the position

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

get_lon(X::SRID) = X.x
get_lat(X::SRID) = X.y
get_alt(X::SRID) = X.z

get_east(X::SRID) = X.x
get_north(X::SRID) = X.y
get_up(X::SRID) = X.z












