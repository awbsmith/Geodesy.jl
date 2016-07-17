
# a handler for Proj4
immutable Proj4Handler   <: AbstractPackageHandler; end     # for point handled by the Proj4 package



"""
Points with a full coordinate reference system as defined by an SRID identifier

Its up to the user to determine the what the x / y / z fields actually represent; which is governed by the element order in Proj4

For a quick reference:
    lat long style CRS's,  x -> lon, y -> lat (or get_lat() and get_lon())
    utm style CRS's,  x -> false east, y -> false north, z -> up (or get_east() get_north() get_up())
"""
immutable CRS{T <: AbstractSRID} <: WorldPosition
    x::Float64
    y::Float64
    z::Float64
end

# useful shortcuts
typealias CRS_NULL CRS{UnknownSRID}

default_params{T <: CRS}(::Type{T}) = (UnknownSRID,)


# trait style functions
has_srid{T <: CRS}(::Type{T}) = Val{true}
get_handler{T <: CRS}(::Type{T}) = Proj4Handler




######################################################################################
### Experimental, compound coordinate reference systems (CCRS)                     ###
######################################################################################

"""
Abstract type for compound coordinate reference system (i.e. height is not ellipsoidal)
"""
abstract AbstractCCRS{T, U} <: WorldPosition

# use the SRID style because we need to Proj4 to handle the Geoid anyway
"""
Compound coordinate reference system where the height is relative to a geoid
"""
immutable CCRS_Geoid{T <: AbstractSRID, U <: AbstractGeoid} <: AbstractCCRS{T, U}
    x::Float64
    y::Float64
    z::Float64
end

default_params{T <: CCRS_Geoid}(::Type{T}) = (UnknownSRID,
                                              UnknownGeoid)

# trait style functions
has_srid{T <: CCRS_Geoid}(::Type{T}) = Val{true}
has_geoid{T <: CCRS_Geoid}(::Type{T}) = Val{true}
get_handler{T <: CCRS_Geoid}(::Type{T}) = Proj4Handler


# Use Geodesy to initialize a standard set of methods
Proj4_types = [CRS, CCRS_Geoid]
for (i, t) in enumerate(Proj4_types)
     Geodesy.eval(build_methods(t, fieldnames(t), Geodesy.GeodesyTypes))
end

