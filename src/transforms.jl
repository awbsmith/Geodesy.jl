#
# Generic version of geotransform, fall through to including the handlers as the final two inputs
#
"""
  geotransform - function to transfom points between coordinate systems

  Usage:

    1) geotransform(OutputType, X)

        Transform geodesy point X to the point type defined by OutputType.  If the output type has an unknown coordinate system / datum,
        it will be assigne the same one as the input type

        Examples:
            geotransform(ECEF, LLA(0,0,0,WGS84))

    2) geotransform(OutputType, OutputCRS, X, InputCRS)

        Transform geodesy point X to the coordinate system and datum defined by OutputCRS.  Input CRS is taken from InputCRS and not from X, even if its available in X

        Examples:
            geotransform(ECEF, CRS(ECEF_CS, WGS84), LLA(0,0,0), CRS(LLA_CS, WGS84))


    3) geotransform(OutputType, X, handler)

        Like definition 1, however the transformation will be performed by the package specified by handler

        Examples:
            geotransform(ECEF, LLA(0,0,0,WGS84), GeodesyHandler)  # use this package to perform the transform

    4) geotransform(OutputType, OutputCRS, X, InputCRS, handler)

        Like definition 3, however the transformation will be performed by the package specified by handler

        Examples:
            geotransform(ECEF, CRS(ECEF_CS, WGS84), LLA(0,0,0), CRS(LLA_CS, WGS84), GeodesyHandler)  # use this package to perform the transform



"""

#
# send everything to the full version of the transform
#
@inline geotransform{oT}(::Type{oT}, X) = geotransform(oT, get_crs(oT), X, get_crs(X))
@inline geotransform{oT, Handler}(::Type{oT}, X, ::Type{Handler}) = geotransform(oT, get_crs(oT), X, get_crs(X), Handler)
@inline geotransform{oT, CRS_out, CRS_in}(::Type{oT}, crs_out::CRS_out, X, crs_in::CRS_in) = geotransform(oT, crs_out, X, crs_in, get_handler(CRS_out, CRS_in))

#
# Overloads for unhandled conversions
#
geotransform{oT, CRS_out, CRS_in, Handler}(::Type{oT}, crs_out::CRS_out, X, crs_in::CRS_in, ::Type{Handler}) = error("The package with handler $(Handler) has no transform to $(crs_out) from $(crs_in) defined ")

# allow promoting datums to CRS's if we know we should
DatumTypes = Union{AbstractDatum, AbstractPosition} # allow positions to be used as datums (e.g. for ENU)

# both of the inbuilt handlers split the data the same way
InbuiltHandlers = Union{CoordinateTransformsHandler, GeodesyHandler}

@inline geotransform{oT, Handler <: InbuiltHandlers}(::Type{oT}, datum_out::DatumTypes, X, datum_in::DatumTypes, ::Type{Handler}) =
              geotransform(oT, CRS(get_coord_system(oT), datum_out), X, CRS(get_coord_system(X), datum_in), Handler)

@inline geotransform{oT, Handler <: InbuiltHandlers}(::Type{oT}, datum::DatumTypes, X, crs_in, ::Type{Handler}) =
              geotransform(oT, CRS(get_coord_system(oT), datum), X, crs_in, Handler)

@inline geotransform{oT, Handler <: InbuiltHandlers}(::Type{oT}, crs_out, X, datum::DatumTypes, ::Type{Handler}) =
              geotransform(oT, crs_out, X, CRS(get_coord_system(X), datum), Handler)

# using the geodesy package, split the CRS's into CS and Datums
@inline geotransform{oT, Handler <: InbuiltHandlers}(::Type{oT}, crs_out, X, crs_in, ::Type{Handler}) =
                geotransform(oT, get_coord_system(crs_out), get_datum(crs_out), X, get_coord_system(crs_in), get_datum(crs_in), Handler)

#
# now fill in the methods
#
@inline geotransform{oT, Handler <: InbuiltHandlers}(::Type{oT}, CS1, datum_out, Xin, CS2, datum_in, ::Type{Handler}) = error("The $(Handler) handler does not define a transform from $(CS1) to $(CS2)")


###############################
### LLA to ECEF coordinates ###
###############################

@inline function geotransform{oT}(::Type{oT}, cs_out::ECEF_CS, datum_out, Xin, cs_in::Union{LL_CS, LLA_CS}, datum_in, ::Type{GeodesyHandler})
    elliptic_datum = combine_datums(datum_out, datum_in)
    assemble_position(oT, lla_to_ecef(Xin, get_parameters(datum_ellipse(elliptic_datum))), CRS(cs_out, elliptic_datum))
end

function lla_to_ecef(lla, ellipse)
     ϕdeg, λdeg, h = get_lat(lla), get_lon(lla), get_alt(lla)

    sinϕ, cosϕ = sind(ϕdeg), cosd(ϕdeg)
    sinλ, cosλ = sind(λdeg), cosd(λdeg)

    N = ellipse.a / sqrt(1 - ellipse.e² * sinϕ^2)  # Radius of curvature (meters)

    x = (N + h) * cosϕ * cosλ
    y = (N + h) * cosϕ * sinλ
    z = (N * (1 - ellipse.e²) + h) * sinϕ

    return (x, y, z)
end


###############################
### ECEF to LLA coordinates ###
###############################


@inline function geotransform{oT}(::Type{oT}, cs_out::LLA_CS, datum_out, Xin, cs_in::ECEF_CS, datum_in, ::Type{GeodesyHandler})
    elliptic_datum = combine_datums(datum_out, datum_in)
    ellipsoid = get_parameters(datum_ellipse(elliptic_datum)) # make sure its the Ellipsoid type
    assemble_position(oT, ecef_to_lla(Xin, ellipsoid), CRS(cs_out, elliptic_datum))
end

@inline function geotransform{oT}(::Type{oT}, cs_out::LL_CS, datum_out, Xin, cs_in::ECEF_CS, datum_in, ::Type{GeodesyHandler})
    elliptic_datum = combine_datums(datum_out, datum_in)
    ellipsoid = get_parameters(datum_ellipse(elliptic_datum)) # make sure its the Ellipsoid type
    lla = ecef_to_lla(Xin, ellipsoid)
    @assert (lla[3] <= alt_tolerance()) "Transforming to LL coordinate sytem but the resultant altitude was not zero"
    assemble_position(oT, (lla[1], lla[2]), CRS(cs_out, elliptic_datum))
end

function ecef_to_lla(ecef, ellipse)

    x, y, z = get_x(ecef), get_y(ecef), get_z(ecef)

    p = hypot(x, y)
    θ = atan2(z*ellipse.a, p*ellipse.b)
    λ = atan2(y, x)
     ϕ = atan2(z + ellipse.e′² * ellipse.b * sin(θ)^3, p - ellipse.e²*ellipse.a*cos(θ)^3)

    N = ellipse.a / sqrt(1 - ellipse.e² * sin(ϕ)^2)  # Radius of curvature (meters)
    h = p / cos(ϕ) - N

    return (rad2deg(λ), rad2deg(ϕ), h)
end


# TESTING - compare accuracy of this method to ecef_to_lla using Geographic lib as a ground truth
function ecef_to_lla_test(ecef, ellipse)

    # termination tolerances for convergence
    hTol = 1e-6
    latTol = hTol / ellipse.a

    # go back into Geodetic via iterative method
    X, Y, Z = get_x(ecef), get_y(ecef), get_z(ecef)
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

    return (rad2deg(lon), rad2deg(lat), h)

end



###############################
### ECEF to ENU coordinates ###
###############################


@inline function geotransform{oT}(::Type{oT}, cs_out::ENU_CS, datum_out, Xin, cs_in::ECEF_CS, datum_in, ::Type{GeodesyHandler})

    # get the output reference position
    origin = datum_position(datum_out)
    elliptic_datum = combine_datums(datum_in, get_datum(origin))
    ellipsoid = get_parameters(datum_ellipse(elliptic_datum)) # make sure its the Ellipsoid type

    # check the ellipsoids match
    assemble_position(oT, ecef_to_enu(Xin, refpos, ellipsoid), CRS(cs_out, datum_out))

end

# worker function
function ecef_to_enu(ecef, ll_ref, ellipse)

    ϕdeg, λdeg, h =  get_lat(ll_ref), get_lon(ll_ref), get_alt(ll_ref)

    # get the reference point in ecef
    ecef_ref = lla_to_ecef(LL(get_lon(ll_ref), get_lat(ll_ref)), ellipse)  # omit height here when calculating the reference point in ECEF.  Add it back at the end
    Δx = get_x(ecef) - ecef_ref[1]
    Δy = get_y(ecef) - ecef_ref[2]
    Δz = get_z(ecef) - ecef_ref[3]

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

    return (east, north, up - h)  # I think the height vector is the same direction as the ellipsoidal height at the ENU's origin
end


###############################
### ENU to ECEF coordinates ###
###############################

@inline function geotransform{oT}(::Type{oT}, cs_out::ECEF_CS, datum_out, Xin, cs_in::ENU_CS, datum_in, ::Type{GeodesyHandler})

    # get the output reference poistion
    refpos = datum_position(datum_in)
    ell_datum = combine_datums(get_datum(refpos), datum_out)
    ellipsoid = get_parameters(datum_ellipse(elliptic_datum)) # make sure its the Ellipsoid type

    # check the ellipsoids match
    assemble_position(oT, enu_to_ecef(Xin, refpos, ellipsoid), CRS(cs_out, ell_datum))

end

# worker function
function enu_to_ecef(enu, ll_ref, ellipse)

     ϕdeg, λdeg, h =  get_lat(ll_ref), get_lon(ll_ref), get_alt(ll_ref)

    # Reference
    ecef_ref = lla_to_ecef(LL(get_lon(ll_ref), get_lat(ll_ref)), ellipse)  # omit height here when calculating the reference point in ECEF.  Add it back at the end

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

    X = ecef_ref[1] + Δx
    Y = ecef_ref[2] + Δy
    Z = ecef_ref[3] + Δz

    return (X,Y,Z)
end


##############################
### LLA to ENU coordinates ###
##############################

@inline function geotransform{oT}(::Type{oT}, cs_out::ENU_CS, datum_out, Xin, cs_in::Union{LLA_CS, LL_CS}, datum_in, ::Type{GeodesyHandler})

    # get the output reference poistion
    refpos = datum_position(datum_out)
    ellipse = datum_ellipse(combine_datums(datum_in, get_datum(get_crs(refpos))))
    ellipsoid = get_parameters(ellipse) # make sure its the Ellipsoid type

    # check the ellipsoids match
    assemble_position(oT, ecef_to_enu(ECEF(lla_to_ecef(Xin, ellipsoid)), refpos, ellipsoid), CRS(cs_out, datum_out))

end




################################
### ENU to LLA coordinates   ###
################################


@inline function geotransform{oT}(::Type{oT}, cs_out::Union{LLA_CS, LL_CS}, datum_out, Xin, cs_in::ENU_CS, datum_in, ::Type{GeodesyHandler})

    # get the output reference position
    refpos = datum_position(datum_in)
    ell_datum = combine_datums(get_datum(get_crs(refpos)), datum_out)
    ellipse = datum_ellipse(ell_datum)
    ellipsoid = get_parameters(ellipse) # make sure its the Ellipsoid type

    # check the ellipsoids match
    assemble_position(oT, ecef_to_lla(ECEF(enu_to_ecef(Xin, refpos, ellipsoid)), ellipsoid), CRS(cs_out, ell_datum))

end

