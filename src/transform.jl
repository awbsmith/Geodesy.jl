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

####################################
### Proj4 backed conversions
### Define first so we can overload
####################################

# back lat lon alt representations on predefined ellipsoids with Proj4 Projections
macro p4_ellipse_str()
	exprs = Vector{Any}(0)
	for ell_type in subtypes(Ellipse)
		str = @sprintf("+proj=longlat +a=%0.19f +b=%0.19f +no_defs", ellipsoid(ell_type).a, ellipsoid(ell_type).b)
		proj = Proj4.Projection(str)
		lhs = :(Proj4.Projection(::Type{LLA{$ell_type}}))  # N.B. argurment is a LLA{ell_type} type not a LLA type directly (because lat long is on the projection string)
		rhs = :($proj)
		push!(exprs, Expr(:(=), lhs, rhs))
	end
	return esc(Expr(:block, exprs...))
end
@p4_ellipse_str()

# now add a generic version of the above
Proj4.Projection{T}(X::LLA{T}) = Proj4.Projection(@sprintf("+proj=longlat +a=%0.19f +b=%0.19f +no_defs", ellipsoid(T).a, ellipsoid(T).b))

# and a generic for SRID points 
Proj4.Projection{T}(X::SRID{T}) = Proj4.Projection{T}



# SRID -> SRID
function Proj4.transform{T <: SRID_Types, U <: SRID_Types}(::Type{SRID{T}}, X::SRID{U}) 
	Y = Proj4.transform(Proj4.Projection(U), Proj4.Projection(T), Vector(X))
	return SRID{T}(Y[1], Y[2], Y[3])
end
Base.convert{T <: SRID_Types}(::Type{SRID{T}}, X::SRID) = Proj4.transform(SRID{T}, X)

# LLA -> SRID
function Proj4.transform{T <: SRID_Types, U <: LLA}(::Type{SRID{T}}, X::U, degs::Bool=true)
	Y = Proj4.transform(Proj4.Projection(U), Proj4.Projection(T), [X.lon, X.lat, X.alt], !degs)   # proj4 is lon lat ordering
	return SRID{T}(Y[1], Y[2], Y[3])
end
Base.convert{T <: SRID_Types}(::Type{SRID{T}}, X::LLA) = Proj4.transform(SRID{T}, X)

# SRID <- LLA 
function Proj4.transform{T <: Datum, U}(::Type{LLA{T}}, X::SRID{U}, degs::Bool=true)
	Y = Proj4.transform(Proj4.Projection(U), Proj4.Projection(LLA{T}), [X.x, X.y, X.z], !degs)   
	return LLA{T}(Y[2], Y[1], Y[3])  # proj4 is lon lat ordering
end
Base.convert{T <: Datum}(::Type{LLA{T}}, X::SRID) = Proj4.transform(LLA{T}, X)






##############################
### LL to ECEF coordinates ###
##############################

function Base.convert{T <: Union{LL, LLA}}(::Type{ECEF}, ll::T)
    ϕdeg, λdeg, h = ll.lat, ll.lon, T <: LLA ? ll.alt : 0
    d = ellipsoid(T)

    sinϕ, cosϕ = sind(ϕdeg), cosd(ϕdeg)
    sinλ, cosλ = sind(λdeg), cosd(λdeg)

    N = d.a / sqrt(1 - d.e² * sinϕ^2)  # Radius of curvature (meters)

    x = (N + h) * cosϕ * cosλ
    y = (N + h) * cosϕ * sinλ
    z = (N * (1 - d.e²) + h) * sinϕ

    return ECEF(x, y, z)
end

##############################
### ECEF to LL coordinates ###
##############################

function Base.convert{T}(::Type{LLA{T}}, ecef::ECEF)
    x, y, z = ecef.x, ecef.y, ecef.z
    d = ellipsoid(T)

    p = hypot(x, y)
    θ = atan2(z*d.a, p*d.b)
    λ = atan2(y, x)
    ϕ = atan2(z + d.e′² * d.b * sin(θ)^3, p - d.e²*d.a*cos(θ)^3)

    N = d.a / sqrt(1 - d.e² * sin(ϕ)^2)  # Radius of curvature (meters)
    h = p / cos(ϕ) - N

    return LLA{T}(rad2deg(ϕ), rad2deg(λ), h)
end

# TODO:
# more coercion than conversion?
# what would not having this a conversion mean for viability of LL type?
function Base.convert{T}(::Type{LL{T}}, ecef::ECEF)
    x, y, z = ecef.x, ecef.y, ecef.z
    d = ellipsoid(T)

    p = hypot(x, y)
    θ = atan2(z*d.a, p*d.b)
    λ = atan2(y, x)
    ϕ = atan2(z + d.e′² * d.b * sin(θ)^3, p - d.e²*d.a*cos(θ)^3)

    N = d.a / sqrt(1 - d.e² * sin(ϕ)^2)  # Radius of curvature (meters)

    return LL{T}(rad2deg(ϕ), rad2deg(λ))
end

###############################
### ECEF to ENU coordinates ###
###############################

function ENU(ecef::ECEF, ll_ref::Union{LL, LLA})
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

    return ENU(east, north, up)
end

###############################
### ENU to ECEF coordinates ###
###############################

function ECEF(enu::ENU, ll_ref::Union{LL, LLA})
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

    return ECEF(X, Y, Z)
end

#############################
### LL to ENU coordinates ###
#############################

function ENU{T <: Union{LL, LLA}}(ll::T, ll_ref::T)
    ecef = ECEF(ll)
    return ENU(ecef, ll_ref)
end

#############################
### ENU to LL coordinates ###
#############################

function Base.call{T <: LLA}(::Type{T}, enu::ENU, ll_ref::Union{LL, LLA})
    ecef = ECEF(enu, ll_ref)
    return T(ecef)
end

function Base.call{T <: LL}(::Type{T}, enu::ENU, ll_ref::Union{LL, LLA})
    ecef = ECEF(enu, ll_ref)
    return T(ecef)
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
