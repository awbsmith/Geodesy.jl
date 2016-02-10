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


##################################################################
### Function to perform transformations                        ###
##################################################################

# uses a common function for a single point
transform_point{T <:Union{WorldPosition, WorldSurfacePosition, LocalPosition}, U <: Union{WorldPosition, WorldSurfacePosition, LocalPosition}}(::Type{T}, X::U) = convert(T, X)

# generic version - overload with more efficient specific conversions
function transform_point{T <:Union{WorldPosition, WorldSurfacePosition, LocalPosition}, U <: AbstractVector}(::Type{T}, X::U)
	Xout = Vector{T}(length(X))
	for (i, Xp) in enumerate(X)
		Xout[i] = transform_point(T, Xp)
	end
	return Xout
end


####################################
### Proj4 backed conversions
### Define first so we can overload
####################################

# build methods to get a proj 4 projection for LLA, ECEF, SRID
Proj4.Projection{T <: EllipseDatum}(X::LLA{T}) = lla_ellipse_proj(T)
Proj4.Projection{T <: EllipseDatum}(::Type{LLA{T}}) = lla_ellipse_proj(T)

Proj4.Projection{T <: EllipseDatum}(X::ECEF{T}) = ecef_ellipse_proj(T)
Proj4.Projection{T <: EllipseDatum}(::Type{ECEF{T}}) = ecef_ellipse_proj(T)

Proj4.Projection{T}(X::SRID_Pos{T}) = get_projection(Val{T})
Proj4.Projection{T}(::Type{SRID_Pos{T}}) = get_projection(Val{T})



##############################################
### LLA / ECEF / SRID "datum" conversions  ###
### Done by Proj4                          ###
##############################################

Proj4.transform{T}(::Type{LLA{T}}, X::LLA{T}) = X		   # same datum
function Proj4.transform{T, U}(::Type{LLA{T}}, X::LLA{U})
	Y = Proj4.transform(Proj4.Projection(X), Proj4.Projection(LLA{T}), [X.lon, X.lat, X.alt])
	return LLA{T}(Y[1], Y[2], Y[3])
end
Base.convert{T}(::Type{LLA{T}}, X::LLA) = Proj4.transform(LLA{T}, X)

Proj4.transform{T <: Datum}(::Type{ECEF{T}}, X::ECEF{T}) = X  # same datum
function Proj4.transform{T <: Datum, U <: Datum}(::Type{ECEF{T}}, X::ECEF{U})
	Y = Proj4.transform(Proj4.Projection(X), Proj4.Projection(ECEF{T}), [X.x, X.y, X.z])
	return ECEF{T}(Y[1], Y[2], Y[3])
end
Base.convert{T}(::Type{ECEF{T}}, X::ECEF) = Proj4.transform(ECEF{T}, X)

# SRID -> SRID
Proj4.transform{T}(::Type{SRID_Pos{T}}, X::SRID_Pos{T}) = X
function Proj4.transform{T, U}(::Type{SRID_Pos{T}}, X::SRID_Pos{U}) 
	Y = Proj4.transform(get_projection(Val{U}), get_projection(Val{T}), Vector(X))
	return SRID_Pos{T}(Y[1], Y[2], Y[3])
end
Base.convert{T}(::Type{SRID_Pos{T}}, X::SRID_Pos) = Proj4.transform(SRID_Pos{T}, X)


######################################
### SRID to / from LLA coordinates ###
######################################

# LLA -> SRID
function Proj4.transform{T}(::Type{Val{T}}, X::LLA, degs::Bool=true)  			         # this version is more natural to type
	Y = Proj4.transform(Proj4.Projection(X), get_projection(Val{T}), [X.lon, X.lat, X.alt], !degs)   # proj4 is lon lat ordering
	return SRID_Pos{T}(Y[1], Y[2], Y[3])
end
Proj4.transform{T}(::Type{SRID_Pos{T}}, X::LLA, degs::Bool=true) = Proj4.transform(Val{T}, X, degs) # but this version is the same syntax as the reverse of the transform
Base.convert{T}(::Type{SRID_Pos{T}}, X::LLA) = Proj4.transform(SRID_Pos{T}, X)

# SRID <- LLA 
function Proj4.transform{T <: Datum}(::Type{LLA{T}}, X::SRID_Pos, degs::Bool=true)
	Y = Proj4.transform(Proj4.Projection(X), lla_ellipse_proj(Val{T}), [X.x, X.y, X.z], !degs)   
	return LLA{T}(Y[1], Y[2], Y[3])  
end
Base.convert{T <: Datum}(::Type{LLA{T}}, X::SRID_Pos) = Proj4.transform(LLA{T}, X)


#######################################
### SRID to / from ECEF coordinates ###
#######################################


# ECEF -> SRID
function Proj4.transform{T}(::Type{Val{T}}, X::ECEF, degs::Bool=true)            # this version is more natural to type
	Y = Proj4.transform(Proj4.Projection(X), get_projection(Val{T}), [X.x, X.y, X.z], !degs)   # proj4 is lon lat ordering
	return SRID_Pos{T}(Y[1], Y[2], Y[3])
end
Proj4.transform{T}(::Type{SRID_Pos{T}}, X::ECEF, degs::Bool=true) = Proj4.transform(Val{T}, X, degs) # but this version is the same syntax as the reverse of the transform
Base.convert{T}(::Type{SRID_Pos{T}}, X::ECEF) = Proj4.transform(SRID_Pos{T}, X)

# SRID <- LLA 
function Proj4.transform{T <: Datum}(::Type{ECEF{T}}, X::SRID_Pos, degs::Bool=true)
	Y = Proj4.transform(Proj4.Projection(X), ecef_ellipse_proj(Val{T}), [X.x, X.y, X.z], !degs)   
	return LLA{T}(Y[1], Y[2], Y[3])  
end
Base.convert{T <: Datum}(::Type{ECEF{T}}, X::SRID_Pos) = Proj4.transform(LLA{T}, X)



##############################
### LL to ECEF coordinates ###
##############################

function Base.convert{T <: Datum}(::Type{ECEF{T}}, ll::Union{LL{T}, LLA{T}})
    ϕdeg, λdeg, h = ll.lat, ll.lon, typeof(ll) <: LLA ? ll.alt : 0
    d = ellipsoid(T)

    sinϕ, cosϕ = sind(ϕdeg), cosd(ϕdeg)
    sinλ, cosλ = sind(λdeg), cosd(λdeg)

    N = d.a / sqrt(1 - d.e² * sinϕ^2)  # Radius of curvature (meters)

    x = (N + h) * cosϕ * cosλ
    y = (N + h) * cosϕ * sinλ
    z = (N * (1 - d.e²) + h) * sinϕ

    return ECEF{T}(x, y, z)
end
Base.convert{T <: Datum, U <: Datum}(::Type{ECEF{T}}, ll::Union{LL{U}, LLA{U}}) = convert(ECEF{T}, convert(ECEF{U}, ll))
Base.convert{T <: Datum}(::Type{ECEF}, ll::Union{LL{T}, LLA{T}}) = convert(ECEF{T}, ll)



##############################
### ECEF to LL coordinates ###
##############################

function Base.convert{T <: Datum}(::Type{LLA{T}}, ecef::ECEF{T})
    x, y, z = ecef.x, ecef.y, ecef.z
    d = ellipsoid(T)

    p = hypot(x, y)
    θ = atan2(z*d.a, p*d.b)
    λ = atan2(y, x)
    ϕ = atan2(z + d.e′² * d.b * sin(θ)^3, p - d.e²*d.a*cos(θ)^3)

    N = d.a / sqrt(1 - d.e² * sin(ϕ)^2)  # Radius of curvature (meters)
    h = p / cos(ϕ) - N

	lla = LLA{T}(rad2deg(λ), rad2deg(ϕ), h)

	
    return lla
end
Base.convert{T <: Datum, U <: Datum}(::Type{LLA{T}}, ecef::ECEF{U}) = convert(LLA{T}, convert(LLA{U}, ecef))
Base.convert{T <: Datum}(::Type{LLA}, ecef::ECEF{T}) = convert(LLA{T}, ecef)

# TESTING - compare accuracy to this method
function convert_test(ecef::ECEF{WGS84})

	# termination tolerances for convergence
    hTol = 1e-6
    latTol = hTol / eWGS84.a

	# go back into Geodetic via iterative method
    X, Y, Z = ecef.x, ecef.y, ecef.z
    lon = atan2(Y, X)
    R = sqrt(X*X + Y*Y)

    # initialisation for iterative lat and h solution
    h = 0.0
    lat = atan2(Z,(1.0 - eWGS84.e²) * R)
	
	converged = false
    maxIter = 10
    i=1
    while (!converged) && (i <= maxIter)
        Nlat = eWGS84.a/(sqrt(1.0 - eWGS84.e² * (sin(lat))^2))
        hNew = R/cos(lat) - Nlat
        latNew = atan2(Z, R * (1.0 -  eWGS84.e² * Nlat / (Nlat+hNew)))
        converged = abs(latNew - lat) < latTol && abs(h-hNew) < hTol
        lat = latNew
        h = hNew
        i=i+1
    end
	
	return LLA{WGS84_ELLIPSE}(rad2deg(lon), rad2deg(lat), h)

end


# TODO:
# more coercion than conversion?
# what would not having this a conversion mean for viability of LL type?
function Base.convert{T <: Datum}(::Type{LL{T}}, ecef::ECEF{T})
    x, y, z = ecef.x, ecef.y, ecef.z
    d = ellipsoid(T)

    p = hypot(x, y)
    θ = atan2(z*d.a, p*d.b)
    λ = atan2(y, x)
    ϕ = atan2(z + d.e′² * d.b * sin(θ)^3, p - d.e²*d.a*cos(θ)^3)

    N = d.a / sqrt(1 - d.e² * sin(ϕ)^2)  # Radius of curvature (meters)

    return LL{T}(rad2deg(λ), rad2deg(ϕ))
end
Base.convert{T <: Datum}(::Type{LL}, ecef::ECEF{T}) = convert(LL{T}, ecef)
Base.call{T <: LL}(::Type{T}, ecef::ECEF) = convert(LL, ecef)


###################################
### Pull the reference position from ###
###################################

# methods to pull the reference point from the template
LLA_ref{T}(::ENU{LLA{T}}) = T


###############################
### ECEF to ENU coordinates ###
###############################

function ENU{T <: Datum}(ecef::ECEF{T}, ll_ref::LLA{T})
    ϕdeg, λdeg = ll_ref.lat, ll_ref.lon

    ecef_ref = ECEF(ll_ref)
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

    return ENU{ll_ref}(east, north, up)
end

###############################
### ENU to ECEF coordinates ###
###############################

# convert to ECEF when no LLA reference point is included in the template for the ENU input
function ECEF{T <: Datum}(enu::ENU, ll_ref::Union{LL{T}, LLA{T}})
    
	ϕdeg, λdeg = ll_ref.lat, ll_ref.lon

    ecef_ref = ECEF(ll_ref)

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

    return ECEF{T}(X, Y, Z)
end

# convert to ECEF when the LLA refernce position is included in the ENU template
ECEF{T <: LLA}(enu::ENU{T}) = ECEF(enu, T)


###############################################
### Transformation matrices for ENU -> ECEF ###
###############################################

# no LLA reference provided in the ENU type so a LLA ref must be supplied
# N.B. this will ignore the LLA refernce in ENU template
function transform_matrix{T <: Datum}(dest::Type{ECEF{T}}, lla_ref::LLA{T})
	
	ϕdeg, λdeg = ll_ref.lat, ll_ref.lon

	# Compute rotation matrix
    sinλ, cosλ = sind(λdeg), cosd(λdeg)
    sinϕ, cosϕ = sind(ϕdeg), cosd(ϕdeg)

	# Reference
	ecef_ref = ECEF{T}(ll_ref)

	Tmat= @fsa([-sinλ  -cosλ*sinϕ    cosλ*cosϕ    ecef_ref.x;
			   cosλ  -sinλ*sinϕ    sinλ*cosϕ    ecef_ref.y;
			   0.0    cosϕ           sinϕ            ecef_ref.z
			   0.0    0.0        0.0          1.0])
end


#############################
### LLA to ENU coordinates ###
#############################

function ENU{T  <: Datum}(lla::LLA{T}, lla_ref::LLA{T})
	ecef = ECEF{T}(lla)
	ENU(ecef, lla_ref)
end


#############################
### ENU to LLA coordinates ###
#############################

# LLA reference NOT provided in the ENU type so it needs to be another input
function Base.call{T <: Datum}(::Type{LLA{T}}, enu::ENU, ll_ref::LLA{T})
    ecef = ECEF(enu, ll_ref)
    return LLA(ecef)
end

# LLA reference provided in the ENU type 
function Base.call{T <: Datum}(::Type{LLA{T}}, enu::ENU)
	lla = LLA(enu)
	@assert(datum(lla) == T, "Output has a different datum to the input")
end
function Base.call(::Type{LLA}, enu::ENU)
	lla_ref = LLA_ref(enu)
    ecef = ECEF(enu, lla_ref)
    return LLA(ecef)
end



################################
### LL to ENU Bounds objects ###
################################

# there's not an unambiguous conversion, but for now,
# returning the minimum bounds that contain all points contained
# by the input bounds
function ENU{T <: Union{LL, LLA}}(bounds::Bounds{T}, ll_ref::T = center(bounds))

    max_x = max_y = -Inf
    min_x = min_y = Inf

    xs = [bounds.min_x, bounds.max_x]
    ys = [bounds.min_y, bounds.max_y]
    if bounds.min_y < 0.0 < bounds.max_y
        push!(ys, 0.0)
    end
    ref_x = getX(ll_ref)
    if bounds.min_x < ref_x < bounds.max_x ||
       (bounds.min_x > bounds.max_x && !(bounds.min_x >= ref_x >= bounds.max_x))
        push!(xs, ref_x)
    end

    for x_ll in xs, y_ll in ys
        pt = ENU(T(y_ll, x_ll), ll_ref)
        x, y = getX(pt), getY(pt)

        min_x, max_x = min(x, min_x), max(x, max_x)
        min_y, max_y = min(y, min_y), max(y, max_y)
    end

    return Bounds{ENU}(min_y, max_y, min_x, max_x)
end
