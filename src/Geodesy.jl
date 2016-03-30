
# we generate code in this module, so precompile where possible
VERSION >= v"0.4.0-dev+6521" && __precompile__(true)

module Geodesy


using FixedSizeArrays
using Proj4  

import Base.show, Base.call, Base.convert, Base.+, Base.-, Base.isnan, Base.getindex

export

    #
    # Position types supported by Geodesy directly
    #

    # World coordinate systems
    LLA,
    LL,
    ECEF,

    # Local coordinate systems
    ENU,  

    #
    # Position types supported by Proj4
    #
    CRS,
    CCRS_Geoid,    


    #
    # convenience for using points with a set datum
    #
    LLA_WGS84,   # typealias for LLA{WGS84}
    ECEF_WGS84,  # typealias for ECEF{WGS84}
    
    # no datum info
    LL_NULL,
    LLA_NULL,   
    ECEF_NULL,
    CRS_NULL,  
    ENU_NULL,

    # Other types
    # Bounds,  # I dont want to export something named Bounds, maybe rename to Geobounds?


    #
    # export some datums
    #
    AbstractDatum,
    UnknownDatum,
    WGS84,
    GRS80,
    GDA94,  

    # datum functions
    get_datum,
    get_datums,    # get a list of all known datums
    
    
    #
    # transformation functions
    #
    geotransform,           # transform a point
    geotransform_vector,    # transform a vector of points
    geotransform_params,    # parameters for a transformation

    #
    # SRID related
    #
    AbstractSRID,
    UnknownSRID, 
    SRID,
    get_srid,
    

    #
    # Geoid related
    #
    AbstractGeoid,  
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
    get_up,

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


for f in ["datums", "ellipsoids", "srid", "geodesy_types", "known_srids", "utm", "math_funcs", "type_methods", "transform", "external",  
          "bounds", "vicenty", "distance",                                                          # do these belong here?
          "Proj4_Types/proj4_types", "Proj4_Types/projections", "Proj4_Types/proj4_transforms"]     # these definitely don't
    include("$f.jl")
end

end # module Geodesy
