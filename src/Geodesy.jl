module Geodesy

export
    # Points
    ECEF,
    ENU,
    LL,
    LLA,
	SRID,

    # Other types
    Bounds,
    Datum,
    Ellipsoid,

    # Constants
    WGS84_ELLIPSE,

	# Srid related
	SRID_Types,
	EPSG_Types,
	ESRI_Types,
	@srid_str,
	EPSG_WGS84,    # typealias for the WGS84 (GPS) SRID type

    # Methods
    center,
    distance,
    getX,
    getY,
    getZ,
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

for f in ["datum", "point", "bounds", "transform", "vicenty", "distance"]
    include("$f.jl")
end

end # module Geodesy
