
abstract AbstractSRID

# when we don't know
immutable UnknownSRID <: AbstractSRID end

# SRID as a type (so we can multiple disbatch based on it)
# auth is a symbol, code is an integer
immutable SRID{auth, code} <: AbstractSRID end  # good to have it as a subtype of datum? 
show{auth, code}(io::IO, ::Type{SRID{auth, code}}) = print(io, "$(auth)$(code)")

#=
# Ideally an SRID would be this but the below can't be use as a parameter of another type for reasons, while a (Symbol, Int) tupple can be a parameter zzz
immutable SRID <: AbstractSRID
    auth::Symbol
    code::Int
end
=#














    
