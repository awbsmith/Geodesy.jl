module Geodesy

export

	# Abstract super types for point types
	WorldPosition, 
    WorldSurfacePosition, 
    LocalPosition,

    # Position types
    ECEF,
    ENU,
    LL,
    LLA,
	SRID,
	LLA_WGS84,   # typealias for LLA{WGS84_ELLIPSE}
	ECEF_WGS84,  # typealias for ECEF{WGS84_ELLIPSE}

	# convenience for using local points without a refernce point
	LLA_Null,  

    # Other types
    Bounds,
    Datum,
    Ellipsoid,

    # Known ellipsoids
    WGS84,

	# Srid related
	SRID_LLA_WGS84,      # typealias for the/a WGS84 (GPS) lat lon alt SRID code
	SRID_ECEF_WGS84,     # typealias for the/a WGS84 (GPS) ECEF SRID code   
	srid_params,         # get the authority and code from the type

    # Methods
    center,
    distance,
    getX,
    getY,
    getZ,
	getlat,
	getlon,
    inBounds

    #= Unexported / Experimental
    ETRS89
    NAD83
    ED50
    OSGB36
    NAD27

    decimal2dms
    dms2decimal

    haversine_distance

    boundaryPoint
    onBounds
    =#

for f in ["datum", "point", "bounds", "utm", "transform", "vicenty", "distance"]
    include("$f.jl")
end

end # module Geodesy
