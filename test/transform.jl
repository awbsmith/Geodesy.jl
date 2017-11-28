using Geodesy
using Geodesy: NAD27, OSGB36, Bounds, center, distance
using Base.Test
using Compat
using StaticArrays

#TODO: Break this up / serious rewrite
@testset "Testing Transforms" begin


###################################
### Testing fixed relationships ###
###################################

@testset "Testing hard coded points" begin

    lla = LLA_WGS84(-71.0960, 42.3673, 0)
    lla_ref = LLA_WGS84(-71.09183, 42.36299, 0)

    # LLA -> ECEF
    ecef = ECEF_WGS84(1529073.1560519305, -4465040.019013103, 4275835.339260309)
    @xyz_approx_eq ECEF(lla) ecef

    #LLA -> ENU
    enu = ENU(-343.493749083977, 478.764855466788, -0.027242885224325164)
    @xyz_approx_eq_eps ENU(lla, lla_ref) enu 1e-8
    @xyz_approx_eq_eps ENU{lla_ref}(lla) enu 1e-8

    #LLA -> ENU (with height in the reference)
    enu_wh = ENU(-343.493749083977, 478.764855466788, -0.027242885224325164 - 10.0)
    lla_ref_wh = LLA_WGS84(-71.09183, 42.36299, 10.0)
    @xyz_approx_eq_eps ENU(lla, lla_ref_wh) enu_wh 1e-8
    @xyz_approx_eq_eps ENU{lla_ref_wh}(lla) enu_wh 1e-8

    # ECEF -> ENU
    @xyz_approx_eq_eps ENU(ecef, lla_ref) enu 1e-8
    @xyz_approx_eq_eps ENU{lla_ref}(ecef) enu 1e-8

    # ENU -> LLA
    @xyz_approx_eq_eps LLA(enu, lla_ref) lla 1e-8
    @xyz_approx_eq_eps LLA(ENU{lla_ref}(enu[1], enu[2], enu[3])) lla 1e-8

    # Bounds{LLA} -> Bounds{ENU}
    bounds = Bounds{LLA_WGS84}(-71.1, -71.094, 42.365, 42.3695)
    bounds_enu_ref = Bounds{ENU}(-247.1268196136449, 247.1268196138187, -249.9308954374605, 249.9353534128848)
    @type_approx_eq Bounds{ENU}(bounds) bounds_enu_ref

end


###################################
### Check CRS / Proj4           ###
###################################

@testset "Testing Proj4 interface" begin

    # make sure the above are still valid
    lla = LLA_WGS84(-71.0960, 42.3673, 0)
    ecef_ref = ECEF_WGS84(1529073.1560519305, -4465040.019013103, 4275835.339260309)

    # check projections are as expected
    @test Geodesy.proj4_str(SRID(LLA_WGS84)) == "+proj=longlat +datum=WGS84 +no_defs"
    @test Geodesy.proj4_str(SRID(ECEF_WGS84)) == "+proj=geocent +datum=WGS84 +units=m +no_defs"

    # and make sure they behave as anticipated
    lla_crs = CRS{SRID(lla)}(lla[1], lla[2], lla[3])
    ecef_crs = CRS{SRID(ecef_ref)}(ecef_ref[1], ecef_ref[2], ecef_ref[3])

    # Non SRID -> SRID
    @xyz_approx_eq_eps CRS{SRID(ecef_ref)}(lla_crs) ecef_ref 1e-8
    @xyz_approx_eq_eps CRS{SRID(lla)}(ecef_crs) lla 1e-8

    # SRID -> Non SRID
    @xyz_approx_eq_eps ECEF(LLA{WGS84}(lla_crs[1], lla_crs[2], lla_crs[3])) ECEF{WGS84}(lla) 1e-8
    @xyz_approx_eq_eps ECEF{WGS84}(ecef_crs[1], ecef_crs[2], ecef_crs[3]) ecef_ref 1e-8

end


#######################################
### Testing ellipsoid relationships ###
#######################################

@testset "Testing Ellipse Relationships" begin

    ecef = ECEF(5.953150599314804e6, 1.5951418955072558e6, 1.6403589592409942e6)
    ecef2wgs = LLA{WGS84}(ecef)
    ecef2nad = LLA{NAD27}(ecef)
    ecef2osgb = LLA{OSGB36}(ecef)

    @test getX(ecef2wgs) == getX(ecef2nad)
    @test abs(getY(ecef2wgs) - getY(ecef2nad)) > 1e-4
    @test abs(getZ(ecef2wgs) - getZ(ecef2nad)) > 10

    @test getX(ecef2wgs) == getX(ecef2osgb)
    @test abs(getY(ecef2wgs) - getY(ecef2osgb)) > 1e-4
    @test abs(getZ(ecef2wgs) - getZ(ecef2osgb)) > 10

end


#######################################
### Test UTM zone / bound           ###
#######################################

@testset "Testing UTM stuff" begin

    srid = SRID{:EPSG, 32755}     # WGS 84 / UTM zone 55 S (approx Sydney, Australia)
    const UTM55S = CRS{srid}    # type alias for convenience

    # get the cnrs (in LLA)
    bounds = Bounds{LLA_WGS84}(144.0, 150.0, -80.0, 0.0)
    dlon = bounds.max_x - bounds.min_x
    dlat = bounds.max_y - bounds.min_y

    # random points inside
    randLLA() = bounds.min_x + rand() * dlon, bounds.min_y + rand() * dlat, (rand() - .5) * 18000

    for i = 1:1000

        # generate (and store for later)
        lla = LLA_WGS84(randLLA())

        # check the zone and bands are reproted correctly
        (zone, band) = Geodesy.utm_zone(lla)

        # correct band and zone zone?
        @test zone == 55
        @test all(band -> -10 <= band < 0, band)

        # correct SRID? (WGS84 datum only)
        (auth, code) = Geodesy.utm_srid(lla)
        srid_out = SRID{auth, code}
        @test srid_out == srid

    end
end

#############################
### Testing random points ###
#############################

@testset "Testing random points" begin

    # lon lat height
    randLLA() = (rand() - .5) * 360, (rand() - .5) * 178, (rand() - .5) * 18000
    srand(0)

    for _ = 1:10_000

        # a random LLA
        x, y, z = randLLA()

        # TODO: figure out how to make this list and draw randomly from it
        # ellipse = Geodesy.get_datums()

        # get LLA and LL types
        lla = LLA{WGS84}(x, y, z)
        ll = LL{WGS84}(x, y)

        # test the center of the bounds
        lla_bounds = Bounds{LLA{WGS84}}(x - 1, x + 1, y - 1, y + 1)
        ll_bounds = Bounds{LLA{WGS84}}(x - 1, x + 1, y - 1, y + 1)

        # check the center is the center
        @xy_approx_eq center(lla_bounds) lla
        @xy_approx_eq center(ll_bounds) ll

        # geotransform to ecef
        ecefa = ECEF(lla)
        ecef = ECEF(ll)

        # test the round trip accuracy
        @xyz_approx_eq_eps LLA(ecefa) lla 1e-6
        @xy_approx_eq_eps LL(ecefa) ll 1e-6

        # test the LLA height match the Euclidean distance between LL and LLA
        @test distance(ecef, ecefa) ≈ abs(getZ(lla)) atol=1e-8

        # get another LL and LLA point
        y, x, z = randLLA()
        lla2 = LLA{WGS84}(x, y, z)
        ll2 = LL{WGS84}(x, y)

        # null ENU point
        enu000 = ENU(0.0, 0.0, 0.0)

        # ecefa is the same point as lla so...
        @xyz_approx_eq_eps ENU(ecefa, lla) enu000   1e-8
        @xy_approx_eq_eps ENU(ecefa, ll) enu000 1e-8
        @z_approx_eq_eps ENU(ecefa, ll) lla     1e-8


        # get another LL and LLA point to use a reference for ENU transforms
        y, x, z = randLLA()
        lla2 = LLA{WGS84}(x, y, z)
        ll2 = LL{WGS84}(x, y)

        # ECEF version of the point
        ecefa2 = ECEF(lla2)
        ecef2 = ECEF(ll2)

        # test LLA -> ECEF -> ENU vs LLA -> ENU
        enu2 = ENU(ecefa, lla2)
        @xyz_approx_eq enu2 ENU(lla, lla2)
        @xy_approx_eq_eps enu2 ENU(ecefa, ll2) 1e-8
        zdiff = getZ(ENU(ecefa, ll2)) - getZ(enu2)
        @test getZ(lla2) ≈ zdiff atol=1e-8

        # Test the transformation matrix approach as well
        (Rf,tf) = Geodesy.geotransform_params(ENU, ECEF, lla2)
        @xyz_approx_eq_eps ENU(Rf * SVector(ecefa) + tf) enu2 1e-8

        # ECEF => ENU => ECEF w/ little change
        enu2v1 = ENU(ecef2, lla)
        @xyz_approx_eq_eps ECEF(enu2v1, lla) ecef2 1e-8

        # Test the transformation matrix approach as well
        (Rb,tb) = Geodesy.geotransform_params(ECEF, ENU, lla)
        @xyz_approx_eq_eps ECEF(Rb * SVector(enu2v1) + tb) ecef2 1e-8

        # ENU => LL same as ENU => ECEF => LLA
        ecef2v1 = ECEF(enu2v1, lla)
        @xyz_approx_eq LLA(enu2v1, lla) LLA(ecef2v1)
        @xy_approx_eq LL(enu2v1, lla) LL(ecef2v1)


        # test distance functions (IDK)
        dist_ecefa = distance(ecefa, ecefa2)
        dist_enua = distance(ENU(lla, lla), ENU(lla2, lla))
        @test dist_ecefa ≈ dist_enua

        dist_ecef = distance(ecef, ecef2)
        dist_enu = distance(ENU(ll, ll), ENU(ll2, ll))
        @test dist_ecef ≈ dist_enu

    end
end

#############################
### Testing vectorization ###
#############################

@testset "Testing Vectorization" begin

    # Step 1 - Grab some points
    lla_data = [147.80465755709005 -35.35851810277833 277.427 "EPSG"  4326;
                147.80051130896    -35.35195210962327 273.8   "EPSG"  4326;
                147.80010482270845 -35.36183505613432 267.334 "EPSG"  4326;
                147.81089289501872 -35.36325561300782 287.688 "EPSG"  4326]

    # test vectorized version vs non vectorized
    Xin = convert(Vector{LLA{WGS84}}, lla_data[:, 1:3]; row=true) # = convert(Vector{CRS{srid}}, raw_data[:, 1:3]; row=true)
    ecef = geotransform(ECEF{WGS84}, Xin)
    for i = 1:length(Xin)
        @xyz_approx_eq ecef[i] ECEF(Xin[i])
    end

    # test vectorization with Proj4
    auth = Symbol(lla_data[1, end-1])
    code = lla_data[1, end]
    srid = SRID{auth, code}

    Xin = convert(Vector{CRS{srid}}, lla_data[:, 1:3]; row=true) # = convert(Vector{CRS{srid}}, raw_data[:, 1:3]; row=true)
    ecef = geotransform(ECEF{WGS84}, Xin)
    for i = 1:length(Xin)
        @xyz_approx_eq ecef[i] ECEF{WGS84}(Xin[i])
    end
end

end


# create a bunch of LLA points for testing
function LLATestSVector(n::Int=1_000_000)

    # get the cnrs (in LLA)
    bounds = Geodesy.Bounds{LLA_WGS84}(144.0, 150.0, -80.0, 0.0)  # utm zone 55
    dlon = bounds.max_x - bounds.min_x
    dlat = bounds.max_y - bounds.min_y

    lla_vec = Vector{LLA{WGS84}}(n)
    for i = 1:n
        lla_vec[i] = LLA{WGS84}(bounds.min_x + rand() * dlon, bounds.min_y + rand() * dlat, (rand() - .5) * 18000)
    end
    return lla_vec

end


function benchmark()

    lla_vec = LLATestSVector(50_000_000)  # random variations have been about ~100ms for this many points

    # make proj4 do it for comparison (so there's no chance for type inference anywhere)
    function P4_conv(lla_vec)

        mat = zeros(length(lla_vec), 3)
        @inbounds for i = 1:length(lla_vec)
            mat[i,1], mat[i,2], mat[i,3] = lla_vec[i][1], lla_vec[i][2], lla_vec[i][3]
        end
        Proj4.geotransform!(Geodesy.get_projection(LLA{WGS84}), Geodesy.get_projection(ECEF{WGS84}), mat)

        ecef_vec = Vector{ECEF{WGS84}}(length(lla_vec))
        @inbounds for i = 1:length(lla_vec)
            ecef_vec[i] = ECEF{WGS84}(mat[i,1], mat[i,2], mat[i,3])
        end
        return ecef_vec

    end

    # make the raw function do it (so there's no chance for type inference anywhere)
    function raw_conv(lla_vec)

        d = Geodesy.ellipsoid(WGS84)

        ecef_vec = Vector{ECEF{WGS84}}(length(lla_vec))
        @inbounds for i = 1:length(lla_vec)

            ϕdeg, λdeg, h = lla_vec[i].lat, lla_vec[i].lon, lla_vec[i].alt

            sinϕ, cosϕ = sind(ϕdeg), cosd(ϕdeg)
            sinλ, cosλ = sind(λdeg), cosd(λdeg)

            N = d.a / sqrt(1 - d.e² * sinϕ^2)  # Radius of curvature (meters)

            x = (N + h) * cosϕ * cosλ
            y = (N + h) * cosϕ * sinλ
            z = (N * (1 - d.e²) + h) * sinϕ


            ecef_vec[i] = ECEF{WGS84}(x,y,z)
        end
        return ecef_vec

    end


    # compile stuff
    raw_conv(lla_vec[1:1])
    geotransform(ECEF, lla_vec[1:1])
    geotransform(ECEF{WGS84}, lla_vec[1:1])
    P4_conv(lla_vec[1:1])

    # and test
    println("Convert to ECEF with no type info in the transform")
    @time raw_conv(lla_vec)

    println("Convert to ECEF")
    @time geotransform(ECEF, lla_vec)

    println("Convert to ECEF{WGS84}")
    @time geotransform(ECEF{WGS84}, lla_vec)

    println("Convert to ECEF{WGS84} via Proj4 (matrix style)")
    @time P4_conv(lla_vec)

end

