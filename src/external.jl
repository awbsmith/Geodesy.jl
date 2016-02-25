##############################################
### Interactions with other point types    ###
##############################################

# convert to a geodesy type - overload this with user methods
geodify{T <: Geodesy_fam}(X::T) = X
geodify(X) = error("No conversion is defined to convert from type $(typeof(X)) to a geodesy type")

# do nothing for geodesy points
ungeodify{T <: Geodesy_fam}(X::T, Y...) = X
ungeodify(X...) = error("No conversion is defined to convert from a geodesy type to a type $(typeof(X[2]))")

# overload this with user methods
ungeodify(X...) = error("No conversion is defined to convert from a geodesy type to a type $(typeof(X[2]))")

# and set a template transformation
transform{T <: Geodesy_fam}(::Type{T}, X) = transform(T, geodify(X))

# and set a default transformation to the desired output type
transform{T <: Geodesy_fam}(::Type{T}, X...) = ungeodify(transform(T, geodify(X[1])), X...)



