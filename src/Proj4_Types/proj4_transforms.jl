##############################################
### SRID "datum" conversions               ###
### (other point types don't have datums)  ###
### Done by Proj4                          ###
##############################################


# CRS{SRID} -> CRS{SRID}
function geotransform{T,U}(::Type{T}, X::U, ::Type{Proj4Handler}, ::Type{Proj4Handler})
    if (T == U)
        out = T(getX(X), getY(X), getZ(X))
    else
        v = [getX(X), getY(X), getZ(X)]
        Proj4.transform!(get_projection(U), get_projection(T), v, false)
        out = T(v[1], v[2], v[3])
    end
end


# CRS{SRID} -> CRS{SRID}, with srids specified by an extra variable (output SRID, input SRID)
function geotransform{T, U, S1, S2}(::Type{T}, X::U, srids::Tuple{S1, S2}, ::Type{Proj4Handler}, ::Type{Proj4Handler}) 
    v = [getX(X), getY(X), getZ(X)]
    Proj4.transform!(get_projection(srids[2]), get_projection(srids[1]), v)
    out = CRS{UnknownSRID}(v[1], v[2], v[3])
end



# X -> CRS{SRID}
function geotransform{T, U}(::Type{T}, X::U, ::Type{Proj4Handler}, ::Type{GeodesyHandler})
    oT = add_param(T, U)
    if false # get_srid(oT) == get_srid(U)  # not valid since the geoid addition
        out = oT(getX(X), getY(X), getZ(X))  # not actually a geotransform.  Should probably be a convert method?
    else
        v = [getX(X), getY(X), getZ(X)]
        Proj4.transform!(get_projection(U), get_projection(oT), v, false)    
        out = oT(v[1], v[2], v[3])
    end
end

# X -> CRS{SRID}, with srids specified by an extra variable (output SRID, input SRID)
function geotransform{T, U, S1, S2}(::Type{T}, X::U, srids::Tuple{S1, S2}, ::Type{Proj4Handler}, ::Type{GeodesyHandler})
    v = [getX(X), getY(X), getZ(X)]
    Proj4.transform!(get_projection(srids[2]), get_projection(srids[1]), v, false)    
    out = CRS{UnknownSRID}(v[1], v[2], v[3])
end

# CRS{SRID} -> X
function geotransform{T, U}(::Type{T}, X::U, ::Type{GeodesyHandler}, ::Type{Proj4Handler})
    oT = add_param(T)
    if false # get_srid(oT) == get_srid(U)  # not valid since the geoid addition
        out = oT(getX(X), getY(X), getZ(X))  # not actually a geotransform.  Should probably be a convert method?
    else
        v = [getX(X), getY(X), getZ(X)]
        Proj4.transform!(get_projection(U), get_projection(oT), v, false)    
        out = oT(v[1], v[2], v[3])
    end
end

# # CRS{SRID} -> X
function geotransform{T, U, S1, S2}(::Type{T}, X::U, srids::Tuple{S1, S2}, ::Type{GeodesyHandler}, ::Type{Proj4Handler})
    v = [getX(X), getY(X), getZ(X)]
    Proj4.transform!(get_projection(srids[2]), get_projection(srids[1]), v, false)    
    out = CRS{UnknownSRID}(v[1], v[2], v[3])
end



#
# Vectorized forms of the above
#
geotransform_vector{T, U}(::Type{T}, X::Vector{U}, ::Type{Proj4Handler}, ::Type{GeodesyHandler}) = proj4_vectorized(T, X)
geotransform_vector{T, U}(::Type{T}, X::Vector{U}, ::Type{GeodesyHandler}, ::Type{Proj4Handler}) = proj4_vectorized(T, X)
geotransform_vector{T, U}(::Type{T}, X::Vector{U}, ::Type{Proj4Handler}, ::Type{Proj4Handler}) = proj4_vectorized(T, X)


geotransform_vector{T, U, S1, S2}(::Type{T}, X::Vector{U}, srids::Tuple{S1, S2}, ::Type{Proj4Handler}, ::Type{GeodesyHandler}) = proj4_vectorized(T, X, srids)
function geotransform_vector{T, U, S1, S2}(::Type{T}, X::Vector{U}, srids::Tuple{S1, S2}, ::Type{GeodesyHandler}, ::Type{Proj4Handler}) 
    proj4_vectorized(T, X, srids)
end
geotransform_vector{T, U, S1, S2}(::Type{T}, X::Vector{U}, srids::Tuple{S1, S2}, ::Type{Proj4Handler}, ::Type{Proj4Handler}) = proj4_vectorized(T, X, srids)

# vectorizied transform for proj4 (looping calls to proj4 is slow)
function proj4_vectorized{T, U}(::Type{T}, X::Vector{U})

    # convert to a matrix 
    mat = Matrix{Float64}(length(X), 3)
    @inbounds for i = 1:length(X)
        mat[i, 1] = getX(X[i])
        mat[i, 2] = getY(X[i])
        mat[i, 3] = getZ(X[i]) 
    end

    # perform it
    Proj4.transform!(Geodesy.get_projection(U), Geodesy.get_projection(T), mat, false)

    # and assign the output
    X = Vector{T}(length(X))
    @inbounds for i = 1:length(X)
        X[i] = T(mat[i,1], mat[i,2], mat[i,3])
    end
    
    return X
end


# vectorizied transform for proj4 (looping calls to proj4 is slow)
function proj4_vectorized{T, U, S1, S2}(::Type{T}, X::Vector{U}, srids::Tuple{S1, S2})

    # convert to a matrix 
    mat = Matrix{Float64}(length(X), 3)
    @inbounds for i = 1:length(X)
        mat[i, 1] = getX(X[i])
        mat[i, 2] = getY(X[i])
        mat[i, 3] = getZ(X[i]) 
    end

    # perform it
    Proj4.transform!(Geodesy.get_projection(srids[2]), Geodesy.get_projection(srids[1]), mat, false)

    # and assign the output
    X = Vector{T}(length(X))
    @inbounds for i = 1:length(X)
        X[i] = T(mat[i,1], mat[i,2], mat[i,3])
    end
    
    return X
end

