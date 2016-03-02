
abstract AbstractSRID

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


# using a generated function to hopefully only generate one projection per SRID,
# no matter how many transforms we do
@generated function get_projection{ T <: SRID}(::Type{T})
    # add the projection info
    const proj = Proj4.Projection(proj4_str(T))  # const is wishfull
    return :($proj)

end







