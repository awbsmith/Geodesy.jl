##############################################
### Interactions with other point types    ###
##############################################

# convert to a geodesy type - overload this with user methods
"""
Convert to geodesy representations of a custom type

If the input is a Type, a Geodesy type should be returned
If the input is an object, a Geodesy object should be returned
"""
geodify{T <: Geodesy_fam}(X::T) = X
geodify{T <: Geodesy_fam}(::Type{T}) = T
geodify{T}(::Type{T}) = error("Function to retrieve the appropriate geodesy type  for a type $(T) is not defined")
geodify(X) = error("No conversion is defined to convert from type $(typeof(X)) to a geodesy type")

# and set a template transformation
geotransform{T <: Geodesy_fam}(::Type{T}, X) = geotransform(T, geodify(X))

# and set a default transformation to the desired output type
geotransform{T}(::Type{T}, X) = T(geotransform(geodify(T), geodify(X)), X)

# The below was for when I was allowing extra arguments to pass to the output type constructor
# geotransform{T}(::Type{T}, X...) = T(geotransform(geodify(T), geodify(X[1])), X...)





