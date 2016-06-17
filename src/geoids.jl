
# a dictionary to allow the user to set the location of gridshift files
geoid_file_dict = Dict{DataType, UTF8String}()


"""
Function to set / get the geoid directory
"""
function set_geoid_dir(path::AbstractString)
    Geodesy.geodesy_properties.geoid_dir = path
end
get_geoid_dir() = Geodesy.geodesy_properties.geoid_dir


# dont know the Geoid
immutable UnknownGeoid <: GeoidDatum; end
geoid(::UnknownGeoid) = error("Unsure which Geoid to use for an unknown Geoid")


#
# Geoid data types
#
abstract GtxGeoidDatum <: GeoidDatum   # for geoids defined by a gridshift file


immutable GTXGeoid <: GtxGeoidDatum # geoid info in a grid shift file
    gtx_file::ASCIIString
end
geoid_file(::GTXGeoid) = X.gtx_file


#
# Geoid tagged types
#

abstract KnownGtxGeoids <: GtxGeoidDatum

# Australia
immutable AusGeoid09 <: KnownGtxGeoids end
geoid_file_dict[AusGeoid09] = "ausgeoid09.gtx"
geoid_file(::AusGeoid09) = geoid_file_dict[AusGeoid09]

# Ellipsoidal (a.k.a no geoid)
immutable NoGeoid <: KnownGtxGeoids end
geoid_file(::NoGeoid) = ""

