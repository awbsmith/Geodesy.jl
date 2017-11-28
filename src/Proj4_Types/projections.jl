# special exception types
@compat abstract type GeodesyException <: Exception end

# unrecognised SRID authority
@compat struct SridAuthException <: GeodesyException
    msg::String
    code::Int
    handler::String
end

# unrecognised SRID authority
@compat struct SridUnknownException <: GeodesyException
    msg::String
    code::Int
    handler::String
end

# unknown geoid
@compat struct GeoidUnknownException <: GeodesyException
    msg::String
    code::Int
    handler::String
end

# can't find the geoid file
@compat struct GeoidFileException <: GeodesyException
    msg::String
    code::Int                          # file not found (SystemError) seems to put a 2 here
end


# get the Proj4 string for a given SRID
function proj4_str{auth, code}(::Type{SRID{auth, code}})

    dict_sym = Symbol(lowercase(string(auth)))

    local dict
    try # hasfield / isfield?
        dict = getfield(Proj4, dict_sym)
    catch
        throw(SridAuthException("Proj4 does not know the SRID Authority: $(auth).\nPlease overload Geodesy.proj4_str to return a the correct Proj4 string for SRID{$(auth), $(code)}\n" *
                                "Geodesy.proj4_str(::Type{SRID{$(auth), $(code)}}) = <Proj4 projection string>", 1, "Proj4"))
    end

    if !haskey(dict, code)
        throw(SridUnknownException("Proj4 does not know the code $(code) for authority $(auth).\nPlease overload Geodesy.proj4_str to return a the correct Proj4 string for SRID{$(auth), $(code)}\n" *
                                   "Geodesy.proj4_str(::Type{SRID{$(auth), $(code)}}) = <Proj4 projection string>", 1, "Proj4"))
    end
    return dict[code]::String

end


# get the Proj4 string for a given SRID
proj4_str{T <: SRID}(::Type{T}, ::Type{NoGeoid}) = proj4_str(T)

# get the Proj4 string for a given SRID and a given geoid
function proj4_str{T <: SRID, G <: KnownGeoid}(::Type{T}, ::Type{G})

    # get the proj4 string for without the geoid first
    p_str = proj4_str(T)

    # add the geoid info
    geoid = geoid_file(G)
    if (length(dirname(geoid)) == 0)  # no path in the filename
        geoid = joinpath(get_geoid_dir(), geoid)
    end
    if !isfile(geoid)
        throw(GeoidFileException("Can not locate the geoid file for geoid: $(G) (Hint: use Geodesy.set_geoid_dir)", 2))
    end
    p_str = String(p_str * " +geoidgrids=$(geoid)")

end


# general error cases
get_projection{T <: UnknownSRID, U <: UnknownGeoid}(::Type{T}, ::Type{U}) = throw(SridUnknownException("Can't build the Proj4 projection for an unknown SRID / unknown Geoid", -1, "Geodesy"))
get_projection{T <: UnknownSRID}(::Type{T}) = throw(SridUnknownException("Can't build the Proj4 projection for an unknown SRID", -1, "Geodesy"))
get_projection{T <: UnknownSRID, U <: KnownGeoid}(::Type{T}, ::Type{U}) = throw(SridUnknownException("Can't build the Proj4 projection for an unknown SRID", -1, "Geodesy"))
get_projection{T <: SRID, U <: UnknownGeoid}(::Type{T}, ::Type{U}) = throw(GeoidUnknownException("Can't build the Proj4 projection for an unknown Geoid", -1, "Geodesy"))


# using a generated function to hopefully only generate one projection per SRID (hopefully this will produce a static var)
# TODO: check this does what I hope
@generated function get_projection{T <: SRID}(::Type{T})
    # add the projection info
    const proj = Proj4.Projection(proj4_str(T))  # const is wishfull
    return :($proj)
end

# get the projection when there's a geoid involved
@generated function get_projection{T <: SRID, G <: KnownGeoid}(::Type{T}, ::Type{G})
    const proj = Proj4.Projection(proj4_str(T, G))  # const is wishfull
    return :($proj)
end

# when the datum is provided as a tuple
get_projection{T}(X::Tuple{T}) = get_projection(X[1])
get_projection{T, G}(X::Tuple{T,G}) = get_projection(X[1], X[2])
