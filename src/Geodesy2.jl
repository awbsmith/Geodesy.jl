
# we generate code in this module, so precompile where possible
#VERSION >= v"0.4.0-dev+6521" && __precompile__(true)

module Geodesy2

using Compat
using Proj4
using FixedSizeArrays
using CoordinateTransformations

export

    #
    # Position types supported by Geodesy directly
    #

    # Generic - used for everthing
    GeoPosition,

    # World point type
    LLA,
    LL,
    ECEF,

    # Local coordinate systems
    ENU,

    # Unknown Position
    UnknownPosition,

    #
    # export some point types with the CRS filled
    #
    LLA_WGS84,
    LL_WGS84,
    ECEF_WGS84,
    ENU_Unknown,

    # Point methods
    get_crs,

    #
    # CRS stuff
    #
    CRS,
    UnknownCRS,
    get_coord_system,
    get_datum,
    datum_ellipse,    # get the datum ellipse
    datum_position,   # get reference position
    datum_time,       # time when the datum is defined


    #
    # export some datums
    #
    UnknownDatum,
    EllipticDatum,
    UnknownEllipticDatum,
    PositionDatum,
    UnknownPositionDatum,
    WGS84,
    GRS80,
    GDA94,

    # datum functions
    get_datums,    # get a list of all known datums

    # export coordinate systems
    LLA_CS,
    LL_CS,
    ECEF_CS,
    ENU_CS,
    UTM_CS,

    # export some CRS typealias's
    LLA_CRS,
    LL_CRS,
    ECEF_CRS,
    ENU_CRS,

    #
    # Position types supported by Proj4
    #
    #CRS,
    #CCRS_Geoid,

    #
    # transformation functions
    #
    geotransform,           # transform a point
    geotransform_vector,    # transform a vector of points
    geotransform_params,    # parameters for a transformation

    #
    # SRID related
    #
    UnknownSRID,
    SRID,
    SRIDt,
    get_srid,


    #
    # Geoid related
    #
    UnknownGeoid,
    get_geoid,

    # utm related
    utm_zone,

    #
    # accessors
    #
    getX,       # should I really export this?
    getY,
    getZ,
    get_lat,
    get_lon,
    get_alt,
    get_east,
    get_north,
    get_up

    #=

    # Other types
    # Bounds,  # I dont want to export something named Bounds, maybe rename to Geobounds?

    # Methods
    # center,   # dont want to export names this generic
    # distance, # dont want to export names this generic

    # convert to / from geodesy types
    geodify,
    ungeodify,

    decimal2dms,
    dms2decimal
    =#


    #=
    inBounds,
    haversine_distance,
    boundaryPoint,
    onBounds
    =#


for f in ["CS", "ellipses", "datum", "elliptic_datums", "geoids", "CRS", "CRSHandler", "SRID", "position_types", "transforms", "utility", "conversion_constructors", "known_srids"]
    include("$f.jl")
end

# include from other packages (other handlers)
for f in ["coordinate_transforms"]
    include("$f.jl")
end



# make this the very last thing!
include("min_display.jl")

end # module Geodesy
