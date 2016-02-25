module Geodesy

using Proj4 # Proj4 is the workhorse behine the module 

import Base.show, Base.call, Base.convert, Proj4.transform

export

    # Generic position types
	CRS,	
	LLA,
	LL,
	ECEF,
	ENU,

	LLA_WGS84,   # typealias for LLA{WGS84}
	ECEF_WGS84,  # typealias for ECEF{WGS84}

	# convenience for using points without a reference point
	LL_NULL,
	LLA_NULL,   
	ECEF_NULL,  
	ENU_NULL,

    # Other types
    Bounds,

    # Named pseudodatums
    WGS84,
    GRS80,
	GDA94,  

	# transform function
	transform, 

	# Srid related
	SRID,

	# utm related
	utmzone, 

    # Methods
    center,
    distance,

	# convert to / from geodesy types
	geodify,
	ungeodify,

	# accessors
	getX,
    getY,
    getZ,
	get_lat,
	get_lon,
	get_alt,
	get_east,
	get_north,
	get_up,

	decimal2dms,
    dms2decimal
    
	#=
	inBounds,
    haversine_distance,
    boundaryPoint,
    onBounds
    =#


for f in ["datums", "ellipsoids", "srids", "geodesy_types", "known_srids", "utm", "point_methods", "transform", "bounds", "vicenty", "distance", "external"]
    include("$f.jl")
end

end # module Geodesy
