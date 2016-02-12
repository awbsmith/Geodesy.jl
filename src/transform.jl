


### Allow stripping the reference out of the templates         ###
convert{T}(::Type{LLA_NULL}, lla::Type{LLA{T}}) = LLA_NULL(lla.lat, lla.lon, lla.alt)		#
convert{T}(::Type{LL_NULL}, ll::Type{LL{T}}) = LL_NULL(lla.lat, lla.lon)		#
convert{T}(::Type{ECEF_NULL}, lla::Type{ECEF{T}}) = ECEF_NULL(lla.lat, lla.lon, lla.alt)		#  
convert{T}(::Type{ENU_NULL}, enu::Type{ENU{T}}) = ENU_NULL(enu.east, enu.north, enu.up)		# allow stripping the reference position out of the template




# more coercion than conversion?
Base.convert{T}(::Type{LL{T}}, lla::LLA{T}) = LL{T}(lla.lon, lla.lat)
Base.convert{T}(::Type{LL}, lla::LLA{T}) = LL{T}(lla.lon, lla.lat)
Base.convert{T}(::Type{LLA{T}}, ll::LL{T}) = LLA{T}(lla.lon, lla.lat, 0.0)
Base.convert{T}(::Type{LLA}, ll::LL{T}) = LLA{T}(lla.lon, lla.lat, 0.0)




####################################
### Function to perform transformations                                  ###
### Defining this so I can worry about what should be a "convert" later  ###
####################################

# do nothing conversions
transform{T}(::Type{LLA{T}}, X::LLA{T}) = X
transform{T}(::Type{LL{T}}, X::LL{T}) = X		   		   
transform{T}(::Type{ECEF{T}}, X::ECEF{T}) = X
transform{T}(::Type{ENU{T}}, X::ENU{T}) = X


# vector format conversions
function transform{T <: Union{WorldPosition, WorldSurfacePosition, LocalPosition}, U <: AbstractVector}(::Type{T}, X::U)
	Xout = Vector{T}(length(X))
	for (i, Xp) in enumerate(X)
		Xout[i] = transform(T, Xp)
	end
	return Xout
end


#
# add some constructor methods while we're here
#
call{T <: WorldPosition, U <: WorldPosition}(::Type{T}, X::U) = transform(T, X)
call{T <: LocalPosition, U <: WorldPosition}(::Type{T}, X::U) = transform(T, X)
call{T <: WorldPosition, U <: LocalPosition}(::Type{T}, X::U) = transform(T, X)


####################################
### Proj4 backed conversions
### Define first so we can overload
####################################

# build methods to get a proj 4 projection (a coordinate reference system) for stuff we know about
get_projection{T <: Union{WorldPosition, WorldSurfacePosition}}(X::T) = get_projection(get_srid(X))




##############################################
### SRID "datum" conversions               ###
### (other point types don't have datums)  ###
### Done by Proj4                          ###
##############################################

function transform{T <: SRID}(::Type{T}, X::SRID_Pos) 
	Y = Proj4.transform(get_projection(X), get_projection(T), Vector(X))
	return SRID_Pos{T}(Y[1], Y[2], Y[3])
end
transform{T <: SRID}(::Type{SRID_Pos{T}}, X::SRID_Pos) = transform(T, X)

# no convert cases
transform{T}(::Type{T}, X::SRID_Pos{T}) = X  # no convert
transform{T}(::Type{SRID_Pos{T}}, X::SRID_Pos{T}) = X		   		   







######################################
### SRID to / from LLA coordinates ###
######################################

# LLA -> SRID
function transform{T <: SRID}(::Type{T}, X::LLA, degs::Bool=true)  			                                # this version is more natural to type
	Y = Proj4.transform(get_projection(X), get_projection(T), [X.lon, X.lat, X.alt], !degs)                 # proj4 is lon lat ordering
	return SRID_Pos{T}(Y[1], Y[2], Y[3])
end
transform{T <: SRID}(::Type{SRID_Pos{T}}, X::LLA, degs::Bool=true) = Proj4.transform(T, X, degs)      # but this version is the same syntax as the reverse of the transform


# SRID -> LLA 
function transform{T <: LLA}(::Type{T}, X::SRID_Pos, degs::Bool=true)
	Y = Proj4.transform(get_projection(X), get_projection(T), [X.x, X.y, X.z], !degs)   
	return T(Y[1], Y[2], Y[3])  
end




#######################################
### SRID to / from ECEF coordinates ###
#######################################


# ECEF -> SRID
function transform{T <: SRID}(::Type{T},  X::ECEF, degs::Bool=true)                     # this version is more natural to type
	Y = Proj4.transform(get_projection(X), get_projection(T), [X.x, X.y, X.z], !degs)   # proj4 is lon lat ordering
	return SRID_Pos{T}(Y[1], Y[2], Y[3])
end
transform{T <: SRID}(::Type{SRID_Pos{T}}, X::ECEF, degs::Bool=true) = Proj4.transform(T, X, degs)      # but this version is the same syntax as the reverse of the transform


# SRID -> ECEF 
function transform{T <: ECEF}(::Type{T}, X::SRID_Pos, degs::Bool=true)
	Y = Proj4.transform(get_projection(X), get_projection(T), [X.x, X.y, X.z], !degs)   
	return T(Y[1], Y[2], Y[3])  
end





##############################
### LL to ECEF coordinates ###
##############################

function transform{T}(::Type{ECEF{T}}, ll::Union{LL{T}, LLA{T}})
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
transform{T}(::Type{ECEF_NULL}, ll::Union{LL{T}, LLA{T}}) = convert(ECEF_NULL, transform(ECEF{T}, ll))
transform{T}(::Type{ECEF}, ll::Union{LL{T}, LLA{T}}) = transform(ECEF{T}, ll)



##############################
### ECEF to LLA coordinates ###
##############################

function transform{T}(::Type{LLA{T}}, ecef::ECEF{T})
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
transform{T}(::Type{LLA_NULL}, ecef::ECEF{T}) = convert(LLA_NULL, transform(LLA{T}, ecef))
transform(::Type{LLA}, ecef::ECEF) = convert(LLA_NULL, ecef)



# Should these exist?
transform{T}(::Type{LL{T}}, ecef::ECEF) = convert(LL{T}, transform(LLA{T}, ecef))
transform(::Type{LL}, ecef::ECEF) = convert(LL_NULL, transform(LLA_NULL, ecef))





# TESTING - compare accuracy of this method to the normal one 
function transform_test{T}(::Type{LLA{T}}, ecef::ECEF{T})

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

# methods to pull the reference point from the template
LLA_ref{T}(::ENU{T}) = T
ELL_type{T}(::LLA{T}) = T


# make the two parameter forms default to using a point as a reference
transform{T <: LLA}(::Type{ENU{T}}, ecef::ECEF) = ENU{T}(ecef_to_enu(ecef, T)...)                    # user specified, wonderful
transform{T}(::Type{ENU}, ecef::ECEF{T}) = ENU{transform(LLA{T}, ecef)}(0.0, 0.0, 0.0)        # assume the input ecef point is the 0

# return the null reference point variety
transform{T}(::Type{ENU_NULL}, ecef::ECEF{T}, ll_ref::LLA{T}) = ENU_NULL(ecef_to_enu(ecef, ll_ref)...)
transform{T}(::Type{ENU}, ecef::ECEF{T}, ll_ref::LLA{T}) = ENU{ll_ref}(ecef_to_enu(ecef, ll_ref)...)
transform{T <: LLA}(::Type{ENU{T}}, ecef::ECEF, ll_ref::LLA) = error("Don't specify both the reference point in the type and provide it as an argument")


# worker function
function ecef_to_enu{T}(ecef::ECEF{T}, ll_ref::LLA{T})
    ϕdeg, λdeg = ll_ref.lat, ll_ref.lon

    ecef_ref = transform(ECEF{T}, ll_ref)
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



###############################
### ENU to ECEF coordinates ###
###############################

# convert to ECEF when the LLA reference position is included in the ENU template
transform{T, U <: LLA}(::Type{ECEF{T}}, enu::ENU{U}) = transform(ECEF{T}, ECEF{ELL_type{LLA_ref{U}}}(enu_to_ecef(enu, U)...))
transform{T <: LLA}(::Type{ECEF}, enu::ENU{T}) = ECEF{ELL_type{T}}(enu_to_ecef(enu, T)...)


# convert to ECEF when the LLA reference position is not included in the ENU template
transform{T}(::Type{ECEF}, enu::ENU, lla_ref::LLA{T}) = ECEF{T}(enu_to_ecef(enu, lla_ref)...)
transform{T}(::Type{ECEF{T}}, enu::ENU, lla_ref::LLA{T}) = ECEF{T}(enu_to_ecef(enu, lla_ref)...)



# convert to ECEF when no LLA reference point is included in the template for the ENU input
function enu_to_ecef{T}(enu::ENU, ll_ref::Union{LL{T}, LLA{T}})
    
	ϕdeg, λdeg = ll_ref.lat, ll_ref.lon

    ecef_ref = transform(ECEF{T}, ll_ref)

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

    return (X,Y,Z)
end




###############################################
### Transformation matrices for ENU -> ECEF ###
###############################################

# no LLA reference provided in the ENU type so a LLA ref must be supplied
# N.B. this will ignore the LLA refernce in ENU template
function transform_matrix{T}(dest::Type{ECEF{T}}, lla_ref::LLA{T})
	
	ϕdeg, λdeg = ll_ref.lat, ll_ref.lon

	# Compute rotation matrix
    sinλ, cosλ = sind(λdeg), cosd(λdeg)
    sinϕ, cosϕ = sind(ϕdeg), cosd(ϕdeg)

	# Reference
	ecef_ref = ECEF{T}(ll_ref)

	Tmat= @fsa([-sinλ  -cosλ*sinϕ    cosλ*cosϕ    ecef_ref.x;
			     cosλ  -sinλ*sinϕ    sinλ*cosϕ    ecef_ref.y;
			     0.0    cosϕ           sinϕ            ecef_ref.z
			     0.0    0.0         0.0         1.0])

end


##############################
### LLA to ENU coordinates ###
##############################

# make the two parameter forms default to using a point as a reference
transform{T,U}(::Type{ENU{T}}, lla::LLA{U}) = ENU{T}(ecef_to_enu(transform(ECEF{U}, lla), T)...)   # user specified, wonderful
transform{T}(::Type{ENU}, lla::LLA{T}) = ENU{lla}(0.0,0.0,0.0)                                     # assume the input point is the 0

# return the null reference point variety
transform{T}(::Type{ENU_NULL}, lla::LLA{T}, ll_ref::LLA{T}) = ENU_NULL(ecef_to_enu(transform(ECEF{T}, lla), ll_ref)...)
transform{T}(::Type{ENU}, lla::LLA{T}, ll_ref::LLA{T}) = ENU{ll_ref}(ecef_to_enu(transform(ECEF{T}, lla), ll_ref)...)
transform{T <: LLA}(::Type{ENU{T}}, lla::LLA, ll_ref::LLA) = error("Don't specify both the reference point in the type and provide it as an argument")



################################
### ENU to LLA coordinates   ###
################################

# convert to ECEF when the LLA reference position is included in the ENU template
transform{T, U <: LLA}(::Type{LLA{T}}, enu::ENU{U}) = transform(LLA{T}, ECEF{ELL_type{LLA_ref{U}}}(enu_to_ecef(enu, U)...))
transform{T <: LLA}(::Type{LLA}, enu::ENU{T}) = LLA{ELL_type{LLA_ref{T}}}(enu_to_ecef(enu, T)...)


# convert to ECEF when the LLA reference position is not included in the ENU template
transform{T}(::Type{LLA}, enu::ENU, lla_ref::LLA{T}) = transform(LLA{T}, ECEF{T}(enu_to_ecef(enu, lla_ref)...))
transform{T}(::Type{LLA{T}}, enu::ENU, lla_ref::LLA{T}) = transform(LLA{T}, ECEF{T}(enu_to_ecef(enu, lla_ref)...))

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
