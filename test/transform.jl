using Geodesy
using Geodesy: NAD27, OSGB36
using Base.Test
using Compat
using FixedSizeArrays

#############################################
### Decimal <=> Degrees, Minutes, Seconds ###
#############################################

for (decimal, d, m, s) in [(0.013, 0.0, 0.0, 46.8),
                           (-0.013, -0.0, 0.0, 46.8),
                           (-0.263, -0.0, 15.0, 46.8),
                           (-179.51, -179.0, 30.0, 36.0)]
    @test Geodesy.dms2decimal(d, m, s) === decimal
    d2, m2, s2 = Geodesy.decimal2dms(decimal)
    @test d2 === d
    @test m2 === m
    @test_approx_eq s2 s
end

################################################
### Helpers for testing approximate equality ###
################################################

# TODO: Move this to Compat.jl
if VERSION < v"0.4.0-dev+3616"
    fieldnames = names
end

macro type_approx_eq(a, b)
    quote
        @test fieldnames($(esc(a))) == fieldnames($(esc(b)))
        for n in fieldnames($(esc(a)))
            @test_approx_eq $(esc(a)).(n) $(esc(b)).(n)
        end
    end
end

macro xyz_approx_eq(a, b)
    quote
        @test_approx_eq getX($(esc(a))) getX($(esc(b)))
        @test_approx_eq getY($(esc(a))) getY($(esc(b)))
        @test_approx_eq getZ($(esc(a))) getZ($(esc(b)))
    end
end
macro xy_approx_eq(a, b)
    quote
        @test_approx_eq getX($(esc(a))) getX($(esc(b)))
        @test_approx_eq getY($(esc(a))) getY($(esc(b)))
    end
end

macro xyz_approx_eq_eps(a, b, eps)
    quote
        @test_approx_eq_eps getX($(esc(a))) getX($(esc(b))) $(esc(eps))
        @test_approx_eq_eps getY($(esc(a))) getY($(esc(b))) $(esc(eps))
        @test_approx_eq_eps getZ($(esc(a))) getZ($(esc(b))) $(esc(eps))
    end
end
macro xy_approx_eq_eps(a, b, eps)
    quote
        @test_approx_eq_eps getX($(esc(a))) getX($(esc(b))) $(esc(eps))
        @test_approx_eq_eps getY($(esc(a))) getY($(esc(b))) $(esc(eps))
    end
end
macro z_approx_eq_eps(a, b, eps)
    quote
        @test_approx_eq_eps getZ($(esc(a))) getZ($(esc(b))) $(esc(eps))
    end
end

###################################
### Testing fixed relationships ###
###################################

lla = LLA_WGS84(-71.0960, 42.3673, 0)
lla_ref = LLA_WGS84(-71.09183, 42.36299, 0)

# LLA -> ECEF
ecef = ECEF_WGS84(1529073.1560519305, -4465040.019013103, 4275835.339260309)
@xyz_approx_eq ECEF(lla) ecef

#LLA -> ENU
enu = ENU(-343.493749083977, 478.764855466788, -0.027242885224325164)
@xyz_approx_eq_eps ENU(lla, lla_ref) enu 1e-8
@xyz_approx_eq_eps ENU{lla_ref}(lla) enu 1e-8

# ECEF -> ENU
@xyz_approx_eq_eps ENU(ecef, lla_ref) enu 1e-8
@xyz_approx_eq_eps ENU{lla_ref}(ecef) enu 1e-8

# ENU -> LLA
@xyz_approx_eq_eps LLA(enu, lla_ref) lla 1e-8
@xyz_approx_eq_eps LLA(ENU{lla_ref}(enu...)) lla 1e-8

# Bounds{LLA} -> Bounds{ENU}
bounds = Bounds{LLA_WGS84}(-71.1, -71.094, 42.365, 42.3695)
bounds_enu_ref = Bounds{ENU}(-247.1268196136449, 247.1268196138187, -249.9308954374605, 249.9353534128848)
@type_approx_eq Bounds{ENU}(bounds) bounds_enu_ref  


###################################
### Check CRS / Proj4         ###
###################################

# make sure the above are still valid
lla = LLA_WGS84(-71.0960, 42.3673, 0)
ecef_ref = ECEF_WGS84(1529073.1560519305, -4465040.019013103, 4275835.339260309)

# check projections are as expected
@test Geodesy.proj4_str(SRID(LLA_WGS84)) == "+proj=longlat +datum=WGS84 +no_defs"
@test Geodesy.proj4_str(SRID(ECEF_WGS84)) == "+proj=geocent +datum=WGS84 +units=m +no_defs"

# and make sure they behave as anticipated
lla_crs = CRS{SRID(lla)}(lla...)
ecef_crs = CRS{SRID(ecef_ref)}(ecef_ref...)

# Non SRID -> SRID
@xyz_approx_eq_eps CRS{SRID(ecef_ref)}(lla_crs) ecef_ref 1e-8
@xyz_approx_eq_eps CRS{SRID(lla)}(ecef_crs) lla 1e-8

# SRID -> Non SRID
@xyz_approx_eq_eps ECEF(LLA{WGS84}(lla_crs...)) ECEF{WGS84}(lla) 1e-8
@xyz_approx_eq_eps ECEF{WGS84}(ecef_crs...) ecef_ref 1e-8


#######################################
### Testing ellipsoid relationships ###
####################################### 

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


#######################################
### Test UTM zone / bound           ###
#######################################

srid = SRID{:EPSG, 32755}     # WGS 84 / UTM zone 55 S (approx Sydney, Australia)
typealias UTM55S CRS{srid}    # type alias for convenience

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
	@test -10 .<= band .< 0

	# correct SRID? (WGS84 datum only)
	srid_out = Geodesy.utm_srid(lla)
	@test srid_out == srid 

end







###############################################
### Test bounding box methods for LLA       ###
###############################################

shift = 180.0 - center(bounds)[1]   # temporarily make lon centered on 180 (to test wrapping)

# add cnr points
lla_vec = [LLA_WGS84(bounds.min_x + shift, bounds.min_y, 0.0), LLA_WGS84(bounds.max_x + shift, bounds.min_y, 0.0), LLA_WGS84(bounds.max_x + shift, bounds.max_y, 0.0), LLA_WGS84(bounds.min_x + shift, bounds.max_y, 0.0)]

# random points inside
randLLA() = bounds.min_x + shift + rand() * dlon, bounds.min_y + rand() * dlat, (rand() - .5) * 18000

# make random points
len = length(lla_vec)
resize!(lla_vec, len+1000)
for i = len+1:len+1000
	lla_vec[i] = randLLA() 
end

# and bounds them
bounds_est = Bounds(lla_vec)

# remove the shift
bounds_est.min_x = Geodesy.bound_thetad(bounds_est.min_x - shift)
bounds_est.max_x = Geodesy.bound_thetad(bounds_est.max_x - shift)

@type_approx_eq bounds bounds_est






#############################
### Testing random errors ###
#############################


# lon lat height
randLLA() = (rand() - .5) * 360, (rand() - .5) * 178, (rand() - .5) * 18000
srand(0)

for _ = 1:10_000

	# a random LLA
    x, y, z = randLLA()

	# TODO: figure out how to make this list and draw randomly from it
	# ellipse = [subtypes(PsuedoDatum{T}); subtypes(PsuedoDynDatum{T})]

	# get LLA and LL types
    lla = LLA{WGS84}(x, y, z)
    ll = LL{WGS84}(x, y)

	# test the center of the bounds
	lla_bounds = Bounds{LLA{WGS84}}(x - 1, x + 1, y - 1, y + 1)
    ll_bounds = Bounds{LLA{WGS84}}(x - 1, x + 1, y - 1, y + 1)
	
	# check the center is the center
	@xy_approx_eq center(lla_bounds) lla
    @xy_approx_eq center(ll_bounds) ll

	# transform to ecef
    ecefa = ECEF(lla)
    ecef = ECEF(ll)

	# test the round trip accuracy
	@xyz_approx_eq_eps LLA(ecefa) lla 1e-6
    @xy_approx_eq_eps LL(ecefa) ll 1e-6

	# test the LLA height match the Euclidean distance between LL and LLA
    @test_approx_eq_eps distance(ecef, ecefa) abs(getZ(lla)) 1e-8

	# get another LL and LLA point
    y, x, z = randLLA()
    lla2 = LLA{WGS84}(x, y, z)
    ll2 = LL{WGS84}(x, y)

	# null ENU point
    enu000 = ENU(0.0, 0.0, 0.0)   

	# ecefa is the same point as lla so...
    @xyz_approx_eq ENU(ecefa, lla) enu000
    @xy_approx_eq_eps ENU(ecefa, ll) enu000 1e-8
    @z_approx_eq_eps ENU(ecefa, ll) lla 1e-8


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
    @test_approx_eq_eps getZ(lla2) zdiff 1e-8

	# Test the transformation matrix approach as well
	(Rf,tf) = Geodesy.transform_params(ENU, lla2)
	@xyz_approx_eq_eps ENU(Rf * Vec(ecefa) + tf) enu2 1e-8

    # ECEF => ENU => ECEF w/ little change
    enu2v1 = ENU(ecef2, lla)
    @xyz_approx_eq_eps ECEF(enu2v1, lla) ecef2 1e-8

	# Test the transformation matrix approach as well
	(Rb,tb) = Geodesy.transform_params(ECEF, lla)
	@xyz_approx_eq_eps ECEF(Rb * Vec(enu2v1) + tb) ecef2 1e-8

    # ENU => LL same as ENU => ECEF => LLA
    ecef2v1 = ECEF(enu2v1, lla)
    @xyz_approx_eq LLA(enu2v1, lla) LLA(ecef2v1)
    @xy_approx_eq LL(enu2v1, lla) LL(ecef2v1)


	# test distance functions (IDK)
    dist_ecefa = distance(ecefa, ecefa2)
    dist_enua = distance(ENU(lla, lla), ENU(lla2, lla))
    @test_approx_eq dist_ecefa dist_enua

    dist_ecef = distance(ecef, ecef2)
    dist_enu = distance(ENU(ll, ll), ENU(ll2, ll))
    @test_approx_eq dist_ecef dist_enu

end


# create a bunch of LLA points for testing
function LLATestVec(n::Int=1_000_000)

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
	
	lla_vec = LLATestVec(50_000_000)  # random variations have been about ~100ms for this many points

	# make proj4 do it for comparison (so there's no chance for type inference anywhere)
	function P4_conv(lla_vec)
		
		mat = zeros(length(lla_vec), 3)
		@inbounds for i = 1:length(lla_vec)
			mat[i,1], mat[i,2], mat[i,3] = lla_vec[i][1], lla_vec[i][2], lla_vec[i][3]
		end
		Proj4.transform!(Geodesy.get_projection(LLA{WGS84}), Geodesy.get_projection(ECEF{WGS84}), mat)
	
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
	transform(ECEF, lla_vec[1:1])
	transform(ECEF{WGS84}, lla_vec[1:1])
	P4_conv(lla_vec[1:1])

	# and test
	println("Convert to ECEF with no type info in the transform")
	@time raw_conv(lla_vec)

	println("Convert to ECEF")
	@time transform(ECEF, lla_vec)

	println("Convert to ECEF{WGS84}")
	@time transform(ECEF{WGS84}, lla_vec)
	
	println("Convert to ECEF{WGS84} via Proj4 (matrix style)")
	@time P4_conv(lla_vec)

end







