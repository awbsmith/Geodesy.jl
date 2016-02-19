using FixedSizeArrays  # to do maths on points

#  TODO: work out why this module kills the convert methods table

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

# Heights relative to something... (basis for a compond coordinate reference system)
abstract  WorldHeight <: Real





########################
### World locations  ###
########################


# proj4 backs this and is lon lat ordering
"""
Point in Longitude-Latitude-Altitude (LLA) coordinates defined for the specified ellipse

Use LLA_NULL(lon, lat, alt) if you don't want to encode the reference ellipse in the type
"""
immutable LLA{T <: AbstractEllipse} <: WorldPosition
	lon::Float64    
	lat::Float64
    alt::Float64
end

# useful shortcuts
typealias LLA_WGS84 LLA{WGS84}
typealias LLA_NULL LLA{UnknownEllipse}




# Global cartesian coordinate system rotating with the Earth
"""
Cartesian cooridnates for a point on an ellipse

Warning: This is a Cartesian system centered on the ellipse's center and with axis direction specified by the ellipse.  This is not necessarily a TRUE ECEF frame, which would be centered on the Earth's 
         center of mass with axis direction given by the International Reference Pole (IRP) and International Reference Meridian

		 Example:
         the "eAiry" ellipse's center is not the Earth's center of mass, so converting from an eAiry based datum LLA{OSGB36} to ECEF{OSGB36} will not give a true ECEF position.  Use the SRID
		 point type to get a true ECEF position if its required.


Use ECEF_NULL(lon, lat, alt) if you don't want to encode the reference ellipse in the type
"""
immutable ECEF{T <: AbstractEllipse} <: WorldPosition
    x::Float64
    y::Float64
    z::Float64
end

# useful shortcuts
typealias ECEF_WGS84 ECEF{WGS84}
typealias ECEF_NULL ECEF{UnknownEllipse}


### SRID based points (converions will use Proj4)
### N.B, for lat lon style SRIDs,  x -> lon, y -> lat (or getlat() and getlon())
### N.B, for utm style SRIDs,  x -> east, y -> north, z -> up (or geteast() getnorth() getup())
immutable SRID_Pos{T <: SRID} <: WorldPosition
   	x::Float64
    y::Float64
	z::Float64
end

# Don't allow unkown SRIDS
Base.call(::Type{SRID_Pos}, x::Real, y::Real, z::Real) = error("Must specify an SRID")

# convenience for getting the srid from an SRID point
SRID{T}(X::SRID_Pos{T}) = T
SRID{T}(::Type{SRID_Pos{T}}) = T





##########################
### Surface locations  ###
##########################


### Point in Latitude-Longitude (LL) coordinates
# proj4 is lon lat ordering
"""
Point in Longitude-Latitude (LL) coordinates defined for the specified ellipse.  Assume the height above the ellipsoid is 0
"""
immutable LL{T <: AbstractEllipse} <: WorldSurfacePosition
	lon::Float64  # proj 4 is lon lat    
	lat::Float64
	LL(x::Real, y::Real) = new(x, y)  # need to specify a constructor to stop the default constructor overwriting the FixedVectorNoTuple{2, Float64} constructors
end


# useful shortcuts
typealias LL_WGS84 LL{WGS84}
typealias LL_NULL LL{UnknownEllipse}




##########################
### Surface heights  ###
##########################

# TODO something with this
immutable EllipHeight{T <: Ellipsoid} <: WorldHeight
	h::Float64
end

# TODO: custom geoid heights
immutable GeoidHeight{T <: AbstractGeoid} <: WorldHeight
	h::Float64
end



######################################################################################
### Experimental, compound coordinate reference systems (CCRS) 					   ###
### Its not intended to work with these, just transform them to/from other types   ###	
######################################################################################

"""
Abstract type for compound coordinate reference system (i.e. height is not ellipsoidal)
"""
abstract CCRS_Pos{T, U} <: WorldPosition

# use the SRID style because we need to Proj4 to handle the Geoid anyway
immutable Geoidal_SRID{T <: SRID, U <: AbstractGeoid} <: CCRS_Pos{T, U}
	x::Float64
	y::Float64
	z::Float64
end



###############################
### Local Coordinate Frames ###
###############################

# adding a layer of abstraction here to allow for ENU points with no LLA reference included in their template
"""
Unknown reference
"""
immutable UnknownRef <: WorldPosition end  # when we don't want to embed the reference frame in out Local coordinates
show(io::IO, ::Type{UnknownRef}) = print(io, "???")


### Point in East-North-Up (ENU) coordinates
# Local cartesian coordinate system
# Linearized about a reference point
"""
East North Up point.  East and North lie in the reference ellipse's tangent plane at the reference point

Use ENU_NULL(e,n,u) if you don't want to encode the reference point in the type
"""
immutable ENU{T} <: LocalPosition   # T should be either UnknownRef or a LLA position
    east::Float64
    north::Float64
    up::Float64
end

typealias ENU_NULL ENU{UnknownRef}
call(::Type{ENU}, e::Real, n::Real) = ENU_NULL(e,n,0.0)  						    # idk
call(::Type{ENU}, e::Real, n::Real, u::Real) = ENU_NULL(e,n,u)  					# allow default constructuction with no reference position


# to add in the template parameter when its omitted
add_param(::Type{ENU}) = ENU_NULL
add_param{T}(::Type{ENU{T}}) = ENU{T}


#ENU(x, y) = ENU(x, y, 0.0)

# TODO: wanted?
#immutable NED <: LocalPosition
#    north::Float64
#	 east::Float64
#    down::Float64
#end
#NED(x, y) = NED(x, y, 0.0)


# retrieve datums and ellipsoids
ellipsoid{T <: AbstractEllipse}(::Union{LLA{T}, LL{T}, ECEF{T}}) = ellipsoid(T)           # reference ellipsoid for the position










