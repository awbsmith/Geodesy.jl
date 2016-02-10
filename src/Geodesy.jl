module Geodesy

import Base.show, Base.call, Base.convert

export

	# Abstract super types for point types
	StaticWorldPosition, 
	WorldPosition, 
    WorldSurfacePosition, 
    LocalPosition,

    # Generic position types
	LLA,
	StaticLLA,
	ECEF,
	SRID_Pos,
	LL,
	ENU,
	LLA_WGS84,   # typealias for LLA{WGS84}
	ECEF_WGS84,  # typealias for ECEF{WGS84}
	LLA_GDA94,   # typealias for StaticLLA{GDA94}

	# convenience for using local points without a refernce point
	NullPos,  

    # Other types
    Bounds,
    Datum,
    Ellipsoid,

    # Named datums (proably should be in its own module)
	# static
    WGS84,
    GRS80,
	
	# dynamic
	GDA94,  


	# Srid related
	SRID,
	get_srid,              # get the authority and code from the type

    # Methods
    center,
    distance,
    getX,
    getY,
    getZ,
	get_lat,
	get_lon,
	get_alt,
	get_east,
	get_north,
	get_up
    

    
	#=
	inBounds,
    decimal2dms,
    dms2decimal,

    haversine_distance,

    boundaryPoint,
    onBounds
    =#


for f in ["datum", "point", "bounds", "utm", "transform", "vicenty", "distance"]
    include("$f.jl")
end

end # module Geodesy
