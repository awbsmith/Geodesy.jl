
abstract AbstractSRID

# when we don't know
immutable UnknownSRID <: AbstractSRID end

# SRID as a type (so we can multiple disbatch based on it)
# auth is a symbol, code is an integer
immutable SRID{auth, code} <: AbstractSRID end  # good to have it as a subtype of datum? 
show{auth, code}(io::IO, ::Type{SRID{auth, code}}) = print(io, "$(auth)$(code)")

#=
# Ideally an SRID would be this but the below can't be use as a parameter of another type for for reasons, while a (Symbol, Int) tupple can zzz
immutable SRID <: AbstractSRID
    auth::Symbol
    code::Int
end
=#


# get the Proj4 string for a given SRID
function proj4_str{auth, code}(::Type{SRID{auth, code}})

    dict_sym = symbol(lowercase(string(auth)))
    
    local dict
    try # hasfield / isfield? 
        dict = Proj4.(dict_sym)
    catch
        error("Proj4 does not know the SRID Authority: $(auth).\nPlease overload Geodesy.proj4_str to return a the correct Proj4 string for SRID{$(auth), $(code)}\n" * 
              "Geodesy.proj4_str(::Type{SRID{$(auth), $(code)}}) = <Proj4 projection string>")
    end

    if !haskey(dict, code)
        error("Proj4 does not know the code $(code) for authority $(auth).\nPlease overload Geodesy.proj4_str to return a the correct Proj4 string for SRID{$(auth), $(code)}\n" * 
              "Geodesy.proj4_str(::Type{SRID{$(auth), $(code)}}) = <Proj4 projection string>")
    end

    return dict[code]::ASCIIString

end


# get the Proj4 string for a given SRID
function proj4_str{T <: SRID, G <: KnownGeoid}(::Type{T}, ::Type{G})

    # get the proj4 string for without the geoid first
    p_str = proj4_str(T)

    # add the geoid info
    geoid = geoid_file(G)
    if (length(dirname(geoid)) == 0)  # no path in the filename
        geoid = joinpath(get_geoid_dir(), geoid)
    end
    if !isfile(geoid)
        error("Cannot find the geoid file for geoid: $(G)\nExpected location was: $(geoid)\n(Hint: use Geodesy.set_geoid_dir)")
    end
    p_str = ASCIIString(p_str * " +geoidgrids=$(geoid)")

end


# general error cases
get_projection{T <: UnknownSRID, U <: UnknownGeoid}(::Type{T}, ::Type{U}) = error("Can't build the Proj4 projection for an unknown SRID / unknown Geoid")
get_projection{T <: UnknownSRID}(::Type{T}) = error("Can't build the Proj4 projection for an unknown SRID")
get_projection{T <: UnknownSRID, U <: KnownGeoid}(::Type{T}, ::Type{U}) = error("Can't build the Proj4 projection for an unknown SRID")
get_projection{T <: SRID, U <: UnknownGeoid}(::Type{T}, ::Type{U}) = error("Can't build the Proj4 projection for an unknown Geoid")


# using a generated function to hopefully only generate one projection per SRID (hopefully this will produce a static var)
# TODO: check this does what I hope
@generated function get_projection{T <: SRID}(::Type{T})
    # add the projection info
    const proj = Proj4.Projection(proj4_str(T))  # const is wishfull
    return :($proj)
end

@generated function get_projection{T <: SRID, G <: KnownGeoid}(::Type{T}, ::Type{G})
    const proj = Proj4.Projection(proj4_str(T, G))  # const is wishfull
    return :($proj)
end








