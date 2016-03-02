

##############################################
### SRID "datum" conversions               ###
### (other point types don't have datums)  ###
### Done by Proj4                          ###
##############################################


function geotransform{T <: SRID, U <: SRID}(::Type{T}, X::CRS{U}) 
    Y = Proj4.transform(get_projection(U), get_projection(T), Vector(X))
    out = CRS{T}(Y[1], Y[2], Y[3])
end
geotransform{T <: SRID, U <: SRID}(::Type{CRS{T}}, X::CRS{U}) = geotransform(T, X)


# X -> SRID
function geotransform{T <: SRID, U <: Proj4_fam}(::Type{T}, X::U)
    iU = add_param(U)
    if T == SRID(iU)
        out = CRS{T}(X[1], X[2], X[3])  # not actually a geotransform.  Should probably be a convert method?
    else
        Y = Proj4.transform(get_projection(iU), get_projection(T), [X[1], X[2], X[3]], false)    
        out = CRS{T}(Y[1], Y[2], Y[3])
    end
    return out
end
geotransform{T <: SRID}(::Type{CRS{T}}, X::Union{ECEF, LLA}) = Proj4.transform(T, X)      # but this version is the same syntax as the reverse of the geotransform


# SRID -> X
function geotransform{T <: Proj4_fam, U <: SRID}(::Type{T}, X::CRS{U})
    iT = add_param(T)
    if SRID(iT) == U
        out = iT(X[1], X[2], X[3])  # not actually a geotransform.  Should probably be a convert method?
    else
        Y = Proj4.transform(get_projection(U), get_projection(iT), [X[1], X[2], X[3]], false)    
        out = iT(Y[1], Y[2], Y[3])
    end
    return out
end



###############################
### LLA to ECEF coordinates ###
###############################


function lla_to_ecef{T <: ECEF, U <: LL_fam}(::Type{T}, ll::U, d::Ellipsoid)
     ϕdeg, λdeg, h = ll.lat, ll.lon, typeof(ll) <: LLA ? ll.alt : 0.0

    sinϕ, cosϕ = sind(ϕdeg), cosd(ϕdeg)
    sinλ, cosλ = sind(λdeg), cosd(λdeg)

    N = d.a / sqrt(1 - d.e² * sinϕ^2)  # Radius of curvature (meters)

    x = (N + h) * cosϕ * cosλ
    y = (N + h) * cosϕ * sinλ
    z = (N * (1 - d.e²) + h) * sinϕ

    return T(x, y, z)
end

# dont allow datum / ellipse transforms within Geodesy
geotransform{T <: KnownDatum, U <: KnownDatum}(::Type{ECEF{T}}, lla::Union{LLA{U}, LL{U}}) = error("Ellipse / datum transforms should be be done explicitly via CRS point types\nIf you're sure what you're doing is correct you can transform the input position to type LLA_NULL")

# must specify the ellipse somewhere
geotransform{T <: UnknownDatum}(::Type{ECEF{T}}, lla::Union{LLA{T}, LL{T}}) = error("An ellipse must be specified in either the input or the output types")

# other combos
geotransform{T <: KnownDatum}(::Type{ECEF{T}}, lla::Union{LLA{T}, LL{T}}) =  lla_to_ecef(ECEF{T}, lla, ellipsoid(T))
geotransform{T <: KnownDatum}(::Type{ECEF_NULL}, lla::Union{LLA{T}, LL{T}}) = lla_to_ecef(ECEF_NULL, lla, ellipsoid(T))
geotransform{T <: KnownDatum}(::Type{ECEF{T}}, lla::LLA_NULL) = lla_to_ecef(ECEF{T}, lla, ellipsoid(T))
geotransform{T <: KnownDatum}(::Type{ECEF}, lla::Union{LLA{T}, LL{T}}) = lla_to_ecef(ECEF{T}, lla, ellipsoid(T))


###############################
### ECEF to LLA coordinates ###
###############################

function ecef_to_lla{T <: LLA, U <: AbstractDatum}(::Type{T}, ecef::ECEF{U}, d::Ellipsoid)
    x, y, z = ecef.x, ecef.y, ecef.z

    p = hypot(x, y)
    θ = atan2(z*d.a, p*d.b)
    λ = atan2(y, x)
     ϕ = atan2(z + d.e′² * d.b * sin(θ)^3, p - d.e²*d.a*cos(θ)^3)

    N = d.a / sqrt(1 - d.e² * sin(ϕ)^2)  # Radius of curvature (meters)
    h = p / cos(ϕ) - N

    lla = T(rad2deg(λ), rad2deg(ϕ), h)
    
    return lla
end

# dont allow datum / ellipse transforms
geotransform{T <: KnownDatum, U <: KnownDatum}(::Type{LLA{T}}, ecef::ECEF{U}) = error("Ellipse / datum transforms should be be done explicitly via CRS point types\nIf you're sure what you're doing is correct you can transform the input position to type ECEF_NULL")

# must specify the ellipse somewhere
geotransform{T <: UnknownDatum}(::Type{LLA{T}}, ecef::ECEF{T}) = error("An ellipse must be specified in either the input or the output types")

# other combos
geotransform{T <: KnownDatum}(::Type{LLA{T}}, ecef::ECEF{T}) = ecef_to_lla(LLA{T}, ecef, ellipsoid(T))
geotransform{T <: KnownDatum}(::Type{LLA_NULL}, ecef::ECEF{T}) = ecef_to_lla(LLA_NULL, ecef, ellipsoid(T))
geotransform{T <: KnownDatum}(::Type{LLA{T}}, ecef::ECEF_NULL) = ecef_to_lla(LLA{T}, ecef, ellipsoid(T))
geotransform{T <: KnownDatum}(::Type{LLA}, ecef::ECEF{T}) = ecef_to_lla(LLA{T}, ecef, ellipsoid(T))

# Should this exist?
geotransform{T <: AbstractDatum}(::Type{LL{T}}, ecef::ECEF) = geotransform(LL{T}, ecef_to_lla(LLA{T}, ecef, ellipsoid(T)))
geotransform{T <: AbstractDatum}(::Type{LL}, ecef::ECEF{T}) = geotransform(LL{T}, ecef_to_lla(LLA{T}, ecef, ellipsoid(T)))



# TESTING - compare accuracy of this method to ecef_to_lla using Geographic lib as a ground truth
function ecef_to_lla_test{T}(::Type{LLA{T}}, ecef::ECEF{T})

    d = ellipsoid(T)

    # termination tolerances for convergence
    hTol = 1e-6
    latTol = hTol / d.a

    # go back into Geodetic via iterative method
    X, Y, Z = ecef.x, ecef.y, ecef.z
    lon = atan2(Y, X)
    R = sqrt(X*X + Y*Y)

    # initialisation for iterative lat and h solution
    h = 0.0
    lat = atan2(Z,(1.0 - d.e²) * R)
    
    converged = false
    maxIter = 10
    i=1
    while (!converged) && (i <= maxIter)
        Nlat = d.a/(sqrt(1.0 - d.e² * (sin(lat))^2))
        hNew = R/cos(lat) - Nlat
        latNew = atan2(Z, R * (1.0 -  d.e² * Nlat / (Nlat+hNew)))
        converged = abs(latNew - lat) < latTol && abs(h-hNew) < hTol
        lat = latNew
        h = hNew
        i=i+1
    end
    
    return LLA_WGS84(rad2deg(lon), rad2deg(lat), h)

end



###############################
### ECEF to ENU coordinates ###
###############################


# make the two parameter forms default to using a point as a reference
geotransform{T}(::Type{ENU}, ecef::ECEF{T}) = error("must supply the reference point in the ENU type or as an additional input")
geotransform{T <: ENU}(::Type{T}, ecef::ECEF) = add_LL_ref(T)(ecef_to_enu(ecef, LL_ref(T))...)   # user specified, wonderful

# return the null reference point variety
geotransform{T}(::Type{ENU_NULL}, ecef::ECEF{T}, ll_ref::Union{LLA{T}, LL{T}}) = ENU_NULL(ecef_to_enu(ecef, ll_ref)...)
geotransform{T}(::Type{ENU}, ecef::ECEF{T}, ll_ref::Union{LLA{T}, LL{T}}) = ENU_NULL(ecef_to_enu(ecef, ll_ref)...)
geotransform{T <: LLA}(::Type{ENU{T}}, ecef::ECEF, ll_ref::LLA) = error("Don't specify both the reference point in the type and provide it as an argument")


# worker function
function ecef_to_enu{T}(ecef::ECEF{T}, ll_ref::Union{LLA{T}, LL{T}})
    ϕdeg, λdeg = ll_ref.lat, ll_ref.lon

    ecef_ref = geotransform(ECEF{T}, ll_ref)
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

    return (east, north, up)
end

###############################################
### Transformation matrices for ECEF -> ENU ###
###############################################

""" 
Returns the rotation and translation to geotransform an ECEF point to a point
centered on ll_ref in the ENU frame:

(R, t) = geotransform_params(ENU, ll_ref)
X_enu = ENU(R * Vec(X_ecef) + t)  # using FixedSizeArrays
 
"""
function geotransform_params{T}(::Type{ENU}, ll_ref::Union{LLA{T}, LL{T}})
    
    ϕdeg, λdeg = ll_ref.lat, ll_ref.lon

    # Compute rotation matrix
    sinλ, cosλ = sind(λdeg), cosd(λdeg)
    sinϕ, cosϕ = sind(ϕdeg), cosd(ϕdeg)

    # Reference
    ecef_ref = ECEF{T}(ll_ref)

    R = @fsa([-sinλ        cosλ        0.0;
              -cosλ*sinϕ   -sinλ*sinϕ    cosϕ;
               cosλ*cosϕ    sinλ*cosϕ    sinϕ])

    t = -R * Vec(ecef_ref.x, ecef_ref.y, ecef_ref.z)

    return (R, t)
    
end
geotransform_params{T <: ENU}(::Type{T}) = geotransform_params(ENU, LL_ref(T))



###############################
### ENU to ECEF coordinates ###
###############################

# convert to ECEF when the LL reference position is included in the ENU template
geotransform{T, U}(::Type{ECEF{T}}, enu::ENU{U}) = geotransform(ECEF{T}, enu_to_ecef(ECEF{ELL_type(U)}, enu, U))
geotransform{T}(::Type{ECEF}, enu::ENU{T}) = enu_to_ecef(ECEF{ELL_type(T)}, enu, T)


# convert to ECEF when the LL reference position is not included in the ENU template
geotransform{T}(::Type{ECEF}, enu::ENU, ll_ref::Union{LLA{T}, LL{T}}) = enu_to_ecef(ECEF{T}, enu, ll_ref)
geotransform{T}(::Type{ECEF{T}}, enu::ENU, ll_ref::Union{LLA{T}, LL{T}}) = enu_to_ecef(ECEF{T}, enu, ll_ref)
geotransform{T}(::Type{ECEF_NULL}, enu::ENU, ll_ref::Union{LLA{T}, LL{T}}) = enu_to_ecef(ECEF_NULL, enu, ll_ref)



# convert to ECEF when no LLA reference point is included in the template for the ENU input
function enu_to_ecef{T <: ECEF, U}(::Type{T}, enu::ENU, ll_ref::Union{LL{U}, LLA{U}})
    
    ϕdeg, λdeg = ll_ref.lat, ll_ref.lon

    ecef_ref = geotransform(ECEF{U}, ll_ref)

    # Compute rotation matrix
    sinλ, cosλ = sind(λdeg), cosd(λdeg)
    sinϕ, cosϕ = sind(ϕdeg), cosd(ϕdeg)

    # Rᵀ = [-sinλ -cosλ*sinϕ cosλ*cosϕ
    #        cosλ -sinλ*sinϕ sinλ*cosϕ
    #         0.0       cosϕ      sinϕ]
    # Δx, Δy, Δz = Rᵀ * [east, north, up]
    Δx = -sinλ * enu.east + -cosλ*sinϕ * enu.north + cosλ*cosϕ * enu.up
    Δy =  cosλ * enu.east + -sinλ*sinϕ * enu.north + sinλ*cosϕ * enu.up
    Δz =   0.0 * enu.east +       cosϕ * enu.north +      sinϕ * enu.up

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
""" 
Returns the rotation and translation to geotransform a point in the ENU frame
centerd on ll_ref into an ECEF point:

(R, t) = geotransform_params(ECEF, ll_ref)  

X_ecef = ECEF(R * Vec(X_enu) + t)  # using FixedSizeArrays
 
"""
function geotransform_params{T, U <: Union{LLA, LL}}(::Type{ECEF{T}}, ll_ref::U)
    
    ϕdeg, λdeg = ll_ref.lat, ll_ref.lon

    # Compute rotation matrix
    sinλ, cosλ = sind(λdeg), cosd(λdeg)
    sinϕ, cosϕ = sind(ϕdeg), cosd(ϕdeg)

    # Reference
    ecef_ref = ECEF{T}(ll_ref)

    R = @fsa([-sinλ  -cosλ*sinϕ    cosλ*cosϕ;
               cosλ  -sinλ*sinϕ    sinλ*cosϕ;
               0.0    cosϕ            sinϕ])

    t = Vec(ecef_ref.x, ecef_ref.y, ecef_ref.z)

    return (R, t)
    
end
geotransform_params{T}(::Type{ECEF}, ll_ref::Union{LLA{T}, LL{T}}) = geotransform_params(ECEF{T}, ll_ref)


##############################
### LLA to ENU coordinates ###
##############################

# make the two parameter forms default to using a point as a reference
geotransform{T <: Union{LLA, LL}}(::Type{ENU}, lla::T) = error("must supply the reference point in the ENU type or as an additional input")
geotransform{T <: ENU, U <: Union{LLA, LL}}(::Type{T}, lla::U) = add_LL_ref(T)(ecef_to_enu(geotransform(ECEF{ELL_type(T)}, lla), LL_ref(T))...)   # user specified, wonderful


# return the null reference point variety
geotransform{T}(::Type{ENU_NULL}, lla::Union{LLA{T}, LL{T}}, ll_ref::Union{LLA{T}, LL{T}}) = ENU_NULL(ecef_to_enu(geotransform(ECEF{T}, lla), ll_ref)...)
geotransform{T}(::Type{ENU}, lla::Union{LLA{T}, LL{T}}, ll_ref::Union{LLA{T}, LL{T}}) = ENU_NULL(ecef_to_enu(geotransform(ECEF{T}, lla), ll_ref)...)
geotransform{T <: LLA}(::Type{ENU{T}}, lla::LLA, ll_ref::LLA) = error("Don't specify both the reference point in the type and provide it as an argument")





################################
### ENU to LLA coordinates   ###
################################

# convert to LLA when the LLA reference position is included in the ENU template
geotransform{T, U <: ENU}(::Type{LLA{T}}, enu::U) = geotransform(LLA{T}, enu_to_ecef(ECEF{ELL_type(U)}, enu, LL_ref(U)))
geotransform{T <: ENU}(::Type{LLA}, enu::T) = geotransform(LLA{ELL_type(T)}, enu_to_ecef(ECEF{ELL_type(T)}, enu, LL_ref(T)))

# Should this exist?
geotransform{T, U <: ENU}(::Type{LL{T}}, enu::U) = geotransform(LLA{T}, enu_to_ecef(ECEF{ELL_type(U)}, enu, LL_ref(U)))
geotransform{T <: ENU}(::Type{LL}, enu::T) = geotransform(LLA{ELL_type(T)}, enu_to_ecef(ECEF{ELL_type(T)}, enu, LL_ref(T)))


# convert to LLA when the LLA reference position is not included in the ENU template
geotransform{T}(::Type{LLA}, enu::ENU, ll_ref::Union{LLA{T}, LL{T}}) = geotransform(LLA{T}, enu_to_ecef(ECEF{T}, enu, ll_ref))
geotransform{T}(::Type{LLA{T}}, enu::ENU, ll_ref::Union{LLA{T}, LL{T}}) = geotransform(LLA{T}, enu_to_ecef(ECEF{T}, enu, ll_ref))


# Should this exist?
geotransform{T <: AbstractDatum}(::Type{LL{T}}, enu::ENU, ll_ref::Union{LLA{T}, LL{T}}) = geotransform(LL{T}, enu_to_ecef(ECEF{T}, enu, ll_ref))
geotransform{T <: AbstractDatum}(::Type{LL}, enu::ENU, ll_ref::Union{LLA{T}, LL{T}}) = geotransform(LL{T}, enu_to_ecef(ECEF{T}, enu, ll_ref))



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



