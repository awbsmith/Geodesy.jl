
#
# Generic version of geotransform, fall through to including the handlers as the final two inputs
# 
geotransform(X)         = geotransform(X, get_handler(X))
geotransform(X, Y)      = geotransform(X, Y, get_handler(X), get_handler(Y))
geotransform(X, Y, ref) = geotransform(X, Y, ref, get_handler(X), get_handler(Y))  # with a reference point

# error message for when no handler is defined
function geotransform{Handler1}(X, ::Type{Handler1}) 
    tOut = isa(X, DataType) ? X : typeof(X)
    tIn = !isa(Y, DataType) ? typeof(Y) : error("geotransform: expected a point type as the second input, got DataType: $(Y)")
    error("No methods defined to transform from $(tIn) to $(tOut)")
end
function geotransform{Handler1, Handler2}(X, Y, ::Type{Handler1}, ::Type{Handler2}) 
    tOut = isa(X, DataType) ? X : typeof(X)
    tIn = !isa(Y, DataType) ? typeof(Y) : error("geotransform: expected a point type as the second input, got DataType: $(Y)")
    error("No methods defined to transform from $(tIn) to $(tOut)")
end
function geotransform{Handler1, Handler2}(X, Y, ref, ::Type{Handler1}, ::Type{Handler2})
    tOut = isa(X, DataType) ? X : typeof(X)
    tIn = !isa(Y, DataType) ? typeof(Y) : error("geotransform: expected a point type as the second input, got DataType: $(Y)")
    error("No methods defined to transform from $(tIn) to $(tOut) using a reference point")
end


#
# version that return the transform parameters
#
geotransform_params(X)         = geotransform_params(X, get_handler(X))
geotransform_params(X, Y)      = geotransform_params(X, Y, get_handler(X), get_handler(Y))
geotransform_params(X, Y, ref) = geotransform_params(X, Y, ref, get_handler(X), get_handler(Y))  # with a reference point


#
# vectorized version (Proj4 point need special handling)
#
geotransform_vector(X)         = geotransform_params(X, get_handler(X))
geotransform_vector(X, Y)      = geotransform_params(X, Y, get_handler(X), get_handler(Y))
geotransform_vector(X, Y, ref) = geotransform_params(X, Y, ref, get_handler(X), get_handler(Y))  # with a reference point



#
# rules for infering one datum from another
#

# the methods
select_datum(::Type{UnknownDatum}, ::Type{UnknownDatum}) = UnknownDatum
select_datum{T <: KnownDatum}(::Type{T}, ::Type{T}) = T
select_datum{T1 <: KnownDatum, T2 <: KnownDatum}(::Type{T1}, ::Type{T2}) = error("Inputs type and output type have different datums")
select_datum{T <: KnownDatum}(::Type{UnknownDatum}, ::Type{T}) = T
select_datum{T <: KnownDatum}(::Type{T}, ::Type{UnknownDatum}) = T


# from 3 options
select_datum{T,U,V}(::Type{T}, ::Type{U}, ::Type{V}) = select_datum(select_datum(get_datum(T), get_datum(U)), get_datum(V))

# fall through methods
select_datum{T,U}(::Type{T}, ::Type{U}) = select_datum(get_datum(T), get_datum(U))
select_datum{T,U}(::T, ::Type{U})       = select_datum(get_datum(T), get_datum(U))
select_datum{T,U}(::Type{T}, ::U)       = select_datum(get_datum(T), get_datum(U))
select_datum{T,U}(::T, ::U)             = select_datum(get_datum(T), get_datum(U))



##########################################
# Special case, allow LL <-> LLA
##########################################

geotransform{T <: LL, U <: LLA}(::Type{T}, lla::U, ::Type{GeodesyHandler}, ::Type{GeodesyHandler}) = LL{select_datum(add_param(T,U), U)}(lla.lon, lla.lat)
geotransform{T <: LLA, U <: LL}(::Type{T}, ll::U, ::Type{GeodesyHandler}, ::Type{GeodesyHandler})  = LLA{select_datum(add_param(T,U), U)}(lla.lon, lla.lat, 0.0)


###############################
### LLA to ECEF coordinates ###
###############################

# reference ellipse inferred from the in / output types 
function geotransform{T <: ECEF, U <: Union{LL, LLA}}(::Type{T}, X::U, ::Type{GeodesyHandler}, ::Type{GeodesyHandler}) 
    oT = add_param(T, U)  # insert the datum if needed
    datum = select_datum(get_datum(oT), get_datum(U))
    ellipse = ellipsoid(datum)
    lla_to_ecef(oT, X, ellipse)
end

# reference ellipse specified (N.B. output type is assumed to be what the user wants)
function geotransform{T <: ECEF, U <: Union{LL, LLA}}(::Type{T}, X::U, ellipse::Ellipsoid, ::Type{GeodesyHandler}, ::Type{GeodesyHandler}) 
    oT = add_param(T, U)  # insert the datum if needed
    lla_to_ecef(oT, X, ellipse)
end

# worker function
function lla_to_ecef{T}(::Type{T}, ll, ellipse::Ellipsoid)
     ϕdeg, λdeg, h = get_lat(ll), get_lon(ll), get_alt(ll)

    sinϕ, cosϕ = sind(ϕdeg), cosd(ϕdeg)
    sinλ, cosλ = sind(λdeg), cosd(λdeg)

    N = ellipse.a / sqrt(1 - ellipse.e² * sinϕ^2)  # Radius of curvature (meters)

    x = (N + h) * cosϕ * cosλ
    y = (N + h) * cosϕ * sinλ
    z = (N * (1 - ellipse.e²) + h) * sinϕ

    return T(x, y, z)
end



###############################
### ECEF to LLA coordinates ###
###############################

# reference ellipse inferred from the in / output types (always using LL to parameterise atm)
function geotransform{T <: Union{LL, LLA}, U <: ECEF}(::Type{T}, X::U, ::Type{GeodesyHandler}, ::Type{GeodesyHandler}) 
    oT = add_param(T, U)  # insert the datum if needed
    datum = select_datum(get_datum(oT), get_datum(U))
    ellipse = ellipsoid(datum)
    ecef_to_lla(oT, X, ellipse)
end

# reference ellipse specified (N.B. output type is assumed to be what the user wants)
function geotransform{T <: Union{LL, LLA}, U <: ECEF}(::Type{T}, X::U, ellipse::Ellipsoid, ::Type{GeodesyHandler}, ::Type{GeodesyHandler}) 
    oT = add_param(T, U)  # insert the datum if needed
    lla_to_ecef(oT, X, ellipse)
end

# worker function
function ecef_to_lla{T <: Union{LL, LLA}}(::Type{T}, ecef, ellipse::Ellipsoid)
    x, y, z = getX(ecef), getY(ecef), getZ(ecef)

    p = hypot(x, y)
    θ = atan2(z*ellipse.a, p*ellipse.b)
    λ = atan2(y, x)
     ϕ = atan2(z + ellipse.e′² * ellipse.b * sin(θ)^3, p - ellipse.e²*ellipse.a*cos(θ)^3)

    N = ellipse.a / sqrt(1 - ellipse.e² * sin(ϕ)^2)  # Radius of curvature (meters)
    h = p / cos(ϕ) - N

    lla = T(rad2deg(λ), rad2deg(ϕ), h)
    
    return lla
end



# TESTING - compare accuracy of this method to ecef_to_lla using Geographic lib as a ground truth
function ecef_to_lla_test{T <: Union{LL, LLA}}(::Type{T}, ecef, ellipse::Ellipsoid)

    ellipse = ellipsoid(T)

    # termination tolerances for convergence
    hTol = 1e-6
    latTol = hTol / ellipse.a

    # go back into Geodetic via iterative method
    X, Y, Z = getX(ecef), getY(ecef), getZ(ecef) 
    lon = atan2(Y, X)
    R = sqrt(X*X + Y*Y)

    # initialisation for iterative lat and h solution
    h = 0.0
    lat = atan2(Z,(1.0 - ellipse.e²) * R)
    
    converged = false
    maxIter = 10
    i=1
    while (!converged) && (i <= maxIter)
        Nlat = ellipse.a/(sqrt(1.0 - ellipse.e² * (sin(lat))^2))
        hNew = R/cos(lat) - Nlat
        latNew = atan2(Z, R * (1.0 -  ellipse.e² * Nlat / (Nlat+hNew)))
        converged = abs(latNew - lat) < latTol && abs(h-hNew) < hTol
        lat = latNew
        h = hNew
        i=i+1
    end
    
    return T(rad2deg(lon), rad2deg(lat), h)

end



###############################
### ECEF to ENU coordinates ###
###############################


# calling with not enough / too much info on the reference point
geotransform(::Type{ENU}, ecef::ECEF, ::Type{GeodesyHandler}, ::Type{GeodesyHandler}) = error("You must supply the reference point in the ENU type or as an additional input")
geotransform{T <: Union{LL, LLA}}(::Type{ENU{T}}, ecef::ECEF, ll_ref, ::Type{GeodesyHandler}, ::Type{GeodesyHandler}) = error("Don't specify both the reference point in the type and provide it as an argument")

# get the reference location from the input type
geotransform{T <: ENU, U <: ECEF}(::Type{T}, X::U, ::Type{GeodesyHandler}, ::Type{GeodesyHandler}) = ecef_to_enu(T, X, get_refloc(T))

# the reference point is user specified
geotransform(::Type{ENU_NULL}, ecef::ECEF, ll_ref::Union{LLA, LL}, ::Type{GeodesyHandler}, ::Type{GeodesyHandler}) = ecef_to_enu(ENU_NULL, ecef, ll_ref)
geotransform(::Type{ENU}, ecef::ECEF, ll_ref::Union{LLA, LL}, ::Type{GeodesyHandler}, ::Type{GeodesyHandler}) = ecef_to_enu(ENU_NULL, ecef, ll_ref)

# if given an arbitrary reference, convert it to LL
geotransform(::Type{ENU}, ecef::ECEF, ref, ::Type{GeodesyHandler}, ::Type{GeodesyHandler}) = geotransform(T, ecef, LL(ref), GeodesyHandler, GeodesyHandler)
geotransform(::Type{ENU_NULL}, ecef::ECEF, ref, ::Type{GeodesyHandler}, ::Type{GeodesyHandler}) = geotransform(T, ecef, LL(ref), GeodesyHandler, GeodesyHandler)


# worker function
function ecef_to_enu{T, U <: ECEF, V <: Union{LLA, LL}}(::Type{T}, ecef::U, ll_ref::V)
    ϕdeg, λdeg, h =  get_lat(ll_ref), get_lon(ll_ref), get_alt(ll_ref)

    # get the reference point in ecef
    ecef_ref = V <: LL ? geotransform(U, ll_ref) : geotransform(U, LL(ll_ref))  # omit height here when calculating the reference point in ECEF.  Add it back at the end
    Δx = ecef.x - ecef_ref.x
    Δy = ecef.y - ecef_ref.y
    Δz = ecef.z - ecef_ref.z

    # Compute rotation matrix
    sinλ, cosλ = sind(λdeg), cosd(λdeg)
    sinϕ, cosϕ = sind(ϕdeg), cosd(ϕdeg)

    # R = [     -sinλ       cosλ  0.0
    #      -cosλ*sinϕ -sinλ*sinϕ cosϕ
    #       cosλ*cosϕ  sinλ*cosϕ sinϕ]
    #
    # east, north, up = R * [Δx, Δy, Δz]
    east  =      -sinλ * Δx +       cosλ * Δy +  0.0 * Δz
    north = -cosλ*sinϕ * Δx + -sinλ*sinϕ * Δy + cosϕ * Δz
    up    =  cosλ*cosϕ * Δx +  sinλ*cosϕ * Δy + sinϕ * Δz

    return T(east, north, up - h)  # I think the height vector is the same direction as the ellipsoidal height at the ENU's origin
end

###############################################
### Transformation matrices for ECEF -> ENU ###
###############################################

""" 
Returns the rotation and translation (R, t) to perform the geotransformation:

Xout = R * Xin + t  # using FixedSizeArrays
 
"""
geotransform_params{T}(::Type{ENU{T}}, ::Type{GeodesyHandler}) = geotransform_params(ENU, T, GeodesyHandler, GeodesyHandler) 

function geotransform_params{T <: Union{LLA, LL}}(::Type{ENU}, ll_ref::T, ::Type{GeodesyHandler}, ::Type{GeodesyHandler})  # TODO: something, because it looks like its transforming LLA -> ENU
    
    ϕdeg, λdeg, h =  get_lat(ll_ref), get_lon(ll_ref), get_alt(ll_ref)

    # Compute rotation matrix
    sinλ, cosλ = sind(λdeg), cosd(λdeg)
    sinϕ, cosϕ = sind(ϕdeg), cosd(ϕdeg)

    # Reference
    ecef_ref = T <: LL ? geotransform(ECEF, ll_ref) : geotransform(ECEF, LL(ll_ref))  # omit height here when calculating the reference point in ECEF.  Add it back at the end

    R = @fsa([-sinλ        cosλ        0.0;
              -cosλ*sinϕ   -sinλ*sinϕ    cosϕ;
               cosλ*cosϕ    sinλ*cosϕ    sinϕ])

    t = -R * Vec(ecef_ref.x, ecef_ref.y, ecef_ref.z) - @fsa([0.0, 0.0, h])
    return (R, t)
    
end



###############################
### ENU to ECEF coordinates ###
###############################

# convert to ECEF when the reference position is included in the ENU template
function geotransform{T <: ECEF, U}(::Type{T}, enu::ENU{U}, ::Type{GeodesyHandler}, ::Type{GeodesyHandler})   # I cant template U into anything because its a value type :-(
    ref = (isa(U, LL) || isa(U, LLA)) ? U : LL(U)
    enu_to_ecef(add_param(T, U), enu, ref)  
end

# convert to ECEF when the LL reference position is not included in the ENU template
geotransform{T <: ECEF}(::Type{T}, enu::ENU_NULL, ll_ref::Union{LLA, LL}, ::Type{GeodesyHandler}, ::Type{GeodesyHandler}) = enu_to_ecef(add_param(T, typeof(ll_ref)), enu, ll_ref)

# if given an arbitrary reference, convert it to LL
geotransform{T <: ECEF}(::Type{T}, enu::ENU_NULL, ref, ::Type{GeodesyHandler}, ::Type{GeodesyHandler}) = geotransform(T, enu, LL(ref), GeodesyHandler, GeodesyHandler)

# spit an error if provided with a reference point in the template and as a parameter
geotransform{T <: ECEF}(::Type{T}, enu::ENU_NULL, ::Type{GeodesyHandler}, ::Type{GeodesyHandler}) = error("You must supply the reference point in the ENU type or as an additional input")
geotransform{T <: ECEF}(::Type{T}, enu::ENU, ll_ref::Union{LLA, LL}, ::Type{GeodesyHandler}, ::Type{GeodesyHandler}) = error("Dont supply both the reference point in thee template and as a parameter")


# worker function
function enu_to_ecef{T <: ECEF, U <: Union{LL, LLA}}(::Type{T}, enu::ENU, ll_ref::U)
    
    ϕdeg, λdeg, h =  get_lat(ll_ref), get_lon(ll_ref), get_alt(ll_ref)

    # Reference
    ecef_ref = U <: LL ? geotransform(T, ll_ref) : geotransform(T, LL(ll_ref))  # omit height here when calculating the reference point in ECEF.  Add it back later

    # Compute rotation matrix
    sinλ, cosλ = sind(λdeg), cosd(λdeg)
    sinϕ, cosϕ = sind(ϕdeg), cosd(ϕdeg)

    # Rᵀ = [-sinλ -cosλ*sinϕ cosλ*cosϕ
    #        cosλ -sinλ*sinϕ sinλ*cosϕ
    #         0.0       cosϕ      sinϕ]
    # Δx, Δy, Δz = Rᵀ * [east, north, up]
    Δx = -sinλ * enu.east + -cosλ*sinϕ * enu.north + cosλ*cosϕ * (enu.up + h)
    Δy =  cosλ * enu.east + -sinλ*sinϕ * enu.north + sinλ*cosϕ * (enu.up + h)
    Δz =   0.0 * enu.east +       cosϕ * enu.north +      sinϕ * (enu.up + h)

    X = ecef_ref.x + Δx
    Y = ecef_ref.y + Δy
    Z = ecef_ref.z + Δz

    return T(X,Y,Z)
end




###############################################
### Transformation matrices for ENU -> ECEF ###
###############################################

# no LLA reference provided in the ENU type so a LLA ref must be supplied
# N.B. this will ignore the LLA refernce in ENU template
function geotransform_params{T <: ECEF, U <: Union{LLA, LL}}(::Type{T}, ll_ref::U, ::Type{GeodesyHandler}, ::Type{GeodesyHandler}) # TODO: something, because it looks like its transforming LLA -> ECEF
    
    ϕdeg, λdeg, h =  get_lat(ll_ref), get_lon(ll_ref), get_alt(ll_ref)

    # Compute rotation matrix
    sinλ, cosλ = sind(λdeg), cosd(λdeg)
    sinϕ, cosϕ = sind(ϕdeg), cosd(ϕdeg)

    # Reference
    ecef_ref = U <: LL ? geotransform(T, ll_ref) : geotransform(T, LL(ll_ref))  # omit height here when calculating the reference point in ECEF.  Add it back later

    R = @fsa([-sinλ  -cosλ*sinϕ    cosλ*cosϕ;
               cosλ  -sinλ*sinϕ    sinλ*cosϕ;
               0.0    cosϕ            sinϕ])

    t = Vec(ecef_ref.x, ecef_ref.y, ecef_ref.z) + R * @fsa([0.0, 0.0, h])

    return (R, t)
    
end



##############################
### LLA to ENU coordinates ###
##############################

# no convert the point to ECEF then go from there
geotransform{T <: ENU, U <: Union{LLA, LL}}(::Type{T}, lla::U, ::Type{GeodesyHandler}, ::Type{GeodesyHandler}) = 
    geotransform(T, ECEF{select_datum(T, U)}(lla), GeodesyHandler, GeodesyHandler)

# convert the point to ECEF then go from there
geotransform{T <: ENU, U <: Union{LLA, LL}, V}(::Type{T}, lla::U, ref::V, ::Type{GeodesyHandler}, ::Type{GeodesyHandler}) = 
    geotransform(T, ECEF{select_datum(T, U, V)}(lla), ref, GeodesyHandler, GeodesyHandler)


################################
### ENU to LLA coordinates   ###
################################

geotransform{T <: Union{LLA, LL}}(::Type{T}, enu::ENU_NULL, ::Type{GeodesyHandler}, ::Type{GeodesyHandler})  = error("You must supply the reference point in the ENU type or as an additional input")

# no convert the point to ECEF then go from there
function geotransform{T <: Union{LLA, LL}, U}(::Type{T}, enu::ENU{U}, ::Type{GeodesyHandler}, ::Type{GeodesyHandler}) 
    eT = ECEF{select_datum(T, U)}  # find a datum somewhere
    T(geotransform(eT, enu, GeodesyHandler, GeodesyHandler))
end

# convert the point to ECEF then go from there
function geotransform{T <: Union{LLA, LL}, U <: ENU, V}(::Type{T}, enu::U, ref::V, ::Type{GeodesyHandler}, ::Type{GeodesyHandler})
    eT = ECEF{select_datum(T, U, V)} # find a datum somewhere
    T(geotransform(eT, enu, ref, GeodesyHandler, GeodesyHandler))
end



#############################################
### Decimal <=> Degrees, Minutes, Seconds ###
#############################################

function decimal2dms(x::Float64)
    d = trunc(x, 0)
    ms = 60 * abs(x - d)
    m = trunc(ms, 0)
    s = 60 * rem(ms, 1)

    return d, m, s
end

function dms2decimal(d::Float64, m::Float64, s::Float64)
    signbit(d) ? d - m/60 - s/3600 :
                 d + m/60 + s/3600
end



#################################################
### A vectorized version of the transforms
#################################################

# vectorizied transform for geodesy types since we need it for Proj4 types
function geotransform_vector{T, U}(::Type{T}, X::Vector{U}, ::Type{GeodesyHandler}, ::Type{GeodesyHandler}) 

    # make sure the output parameter is filled
    oT =  add_param(T, U)
    Xout = Vector{oT}(length(X))
    @inbounds for i = 1:length(X)
        Xout[i] = geotransform(oT, X[i])
    end
    return Xout
end

# vectorizied transform for geodesy types since we need it for Proj4 types
function geotransform_vector{T, U, V}(::Type{T}, X::Vector{U}, ll_ref::V, ::Type{GeodesyHandler}, ::Type{GeodesyHandler}) 

    # make sure the output parameter is filled
    oT =  add_param(T, U)
    Xout = Vector{oT}(length(X))
    @inbounds for i = 1:length(X)
        Xout[i] = geotransform(oT, X[i], ll_ref)
    end
    return Xout
end





##############################################
### SRID "datum" conversions               ###
### (other point types don't have datums)  ###
### Done by Proj4                          ###
##############################################


# CRS{SRID} -> CRS{SRID}
function geotransform{T,U}(::Type{T}, X::U, ::Type{Proj4Handler}, ::Type{Proj4Handler})
    if (get_srid(T) == get_srid(U))
        out = T(getX(X), getY(X), getZ(X))
    else
        Y = Proj4.transform(get_projection(U), get_projection(T), [getX(X), getY(X), getZ(X)])
        out = T(Y[1], Y[2], Y[3])
    end
end


# CRS{SRID} -> CRS{SRID}, with srids specified by an extra variable (output SRID, input SRID)
function geotransform{T, U, S1, S2}(::Type{T}, X::U, srids::Tuple{S1, S2}, ::Type{Proj4Handler}, ::Type{Proj4Handler}) 
    Y = Proj4.transform(get_projection(srids[2]), get_projection(srids[1]), [getX(X), getY(X), getZ(X)])
    out = CRS{UnknownSRID}(Y[1], Y[2], Y[3])
end



# X -> CRS{SRID}
function geotransform{T, U}(::Type{T}, X::U, ::Type{Proj4Handler}, ::Type{GeodesyHandler})
    oT = add_param(T, U)
    if get_srid(oT) == get_srid(U)
        out = oT(getX(X), getY(X), getZ(X))  # not actually a geotransform.  Should probably be a convert method?
    else
        Y = Proj4.transform(get_projection(U), get_projection(oT), [getX(X), getY(X), getZ(X)], false)    
        out = oT(Y[1], Y[2], Y[3])
    end
end

# X -> CRS{SRID}, with srids specified by an extra variable (output SRID, input SRID)
function geotransform{T, U, S1, S2}(::Type{T}, X::U, srids::Tuple{S1, S2}, ::Type{Proj4Handler}, ::Type{GeodesyHandler})
    Y = Proj4.transform(get_projection(srids[2]), get_projection(srids[1]), [getX(X), getY(X), getZ(X)], false)    
    out = CRS{UnknownSRID}(Y[1], Y[2], Y[3])
end

# CRS{SRID} -> X
function geotransform{T, U}(::Type{T}, X::U, ::Type{GeodesyHandler}, ::Type{Proj4Handler})
    oT = add_param(T)
    if get_srid(oT) == get_srid(U)
        out = oT(getX(X), getY(X), getZ(X))  # not actually a geotransform.  Should probably be a convert method?
    else
        Y = Proj4.transform(get_projection(U), get_projection(oT), [getX(X), getY(X), getZ(X)], false)    
        out = oT(Y[1], Y[2], Y[3])
    end
end

# # CRS{SRID} -> X
function geotransform{T, U, S1, S2}(::Type{T}, X::U, srids::Tuple{S1, S2}, ::Type{GeodesyHandler}, ::Type{Proj4Handler})
    Y = Proj4.transform(get_projection(srids[2]), get_projection(srids[1]), [getX(X), getY(X), getZ(X)], false)    
    out = CRS{UnknownSRID}(Y[1], Y[2], Y[3])
end



#
# Vectorized forms of the above
#
geotransform_vector{T, U}(::Type{T}, X::Vector{U}, ::Type{Proj4Handler}, ::Type{GeodesyHandler}) = proj4_vectorized(T, X)
geotransform_vector{T, U}(::Type{T}, X::Vector{U}, ::Type{GeodesyHandler}, ::Type{Proj4Handler}) = proj4_vectorized(T, X)
geotransform_vector{T, U}(::Type{T}, X::Vector{U}, ::Type{Proj4Handler}, ::Type{Proj4Handler}) = proj4_vectorized(T, X)


geotransform_vector{T, U, S1, S2}(::Type{T}, X::Vector{U}, srids::Tuple{S1, S2}, ::Type{Proj4Handler}, ::Type{GeodesyHandler}) = proj4_vectorized(T, X, srids)
geotransform_vector{T, U, S1, S2}(::Type{T}, X::Vector{U}, srids::Tuple{S1, S2}, ::Type{GeodesyHandler}, ::Type{Proj4Handler}) = proj4_vectorized(T, X, srids)
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
    Proj4.transform!(Geodesy.get_projection(S2), Geodesy.get_projection(S1), mat, false)

    # and assign the output
    X = Vector{T}(length(X))
    @inbounds for i = 1:length(X)
        X[i] = T(mat[i,1], mat[i,2], mat[i,3])
    end
    
    return X
end






