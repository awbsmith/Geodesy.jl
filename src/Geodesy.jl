
# we generate code in this module, so precompile where possible
VERSION >= v"0.4.0-dev+6521" && __precompile__(true)

module Geodesy


using FixedSizeArrays
using Proj4  

import Base.show, Base.call, Base.convert

export

    # Position types that understand datums
    CRS,
    CCRS_Geoid,    

    # Position types that understand ellipsoids
    LLA,
    LL,
    ECEF,

    # Local coordinate frames
    ENU,

    LLA_WGS84,   # typealias for LLA{WGS84}
    ECEF_WGS84,  # typealias for ECEF{WGS84}

    # convenience for using points without a reference point
    LL_NULL,
    LLA_NULL,   
    ECEF_NULL,
    CRS_NULL,  
    ENU_NULL,

    # Other types
    # Bounds,  # I dont want to export something named Bounds, maybe rename to Geobounds?

    # Datum related
    get_datum,
    get_datums,    # get a list of all known datums
    AbstractDatum,
    WGS84,
    GRS80,
    GDA94,  
    UnknownDatum,

    # transform function
    geotransform,
    geotransform_vector,  

    # Srid related
    get_srid,
    AbstractSRID,
    SRID,
    UnknownSRID, 

    # Geoid related
    get_geoid,
    AbstractGeoid,  
    UnknownGeoid, 
    

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
    
    
    

    # utm related
    utm_zone, 

    # Methods
    # center,   # dont want to export names this generic
    # distance, # dont want to export names this generic

    # convert to / from geodesy types
    geodify,
    ungeodify,

    decimal2dms,
    dms2decimal

    #=
    inBounds,
    haversine_distance,
    boundaryPoint,
    onBounds
    =#


for f in ["datums", "ellipsoids", "srids", "geodesy_types", "known_srids", "utm", "type_methods", "transform", "bounds", "vicenty", "distance", "external"]
    include("$f.jl")
end

end # module Geodesy
