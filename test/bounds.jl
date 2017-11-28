using Geodesy
using Base.Test
using Geodesy: Bounds, onBounds, inBounds, boundaryPoint, NAD27, center

@testset "Testing Bounds" begin

    # Test bounding box methods for LLA
    @testset "Testing LLA Bounds" begin

        bounds = Bounds{LLA_WGS84}(-1., 1., -1., 1.)
        dlon = bounds.max_x - bounds.min_x
        dlat = bounds.max_y - bounds.min_y

        shift = 180.0 - center(bounds)[1]   # temporarily make lon centered on 180 (to test wrapping)

        # add cnr points
        lla_vec = [LLA_WGS84(bounds.min_x + shift, bounds.min_y, 0.0), LLA_WGS84(bounds.max_x + shift, bounds.min_y, 0.0), LLA_WGS84(bounds.max_x + shift, bounds.max_y, 0.0), LLA_WGS84(bounds.min_x + shift, bounds.max_y, 0.0)]

        # random points inside
        randLLA() = LLA_WGS84(bounds.min_x + shift + rand() * dlon, bounds.min_y + rand() * dlat, (rand() - .5) * 18000)

        # make random points
        append!(lla_vec, [randLLA() for i in 1:1000])

        # and bounds them
        bounds_est = Bounds(lla_vec)

        # remove the shift
        bounds_est.min_x = Geodesy.bound_thetad(bounds_est.min_x - shift)
        bounds_est.max_x = Geodesy.bound_thetad(bounds_est.max_x - shift)

        @type_approx_eq bounds bounds_est
    end

    # Test inbounds
    @testset "Testing inBounds()" for bounds in (Bounds{ENU}(1.1, 2.2, 3.3, 4.4),
                                                 Bounds{LLA{NAD27}}(1.1, 2.2, 3.3, 4.4)
                                                )

        T = Geodesy.point_type(bounds)
        min_x, min_y, max_x, max_y = bounds.min_x, bounds.min_y, bounds.max_x, bounds.max_y

        @test inBounds(T(min_x, min_y), bounds)
        @test !inBounds(T(min_x - eps(min_x), min_y), bounds)
        @test !inBounds(T(min_x, min_y - eps(min_y)), bounds)

        @test inBounds(T(max_x, max_y), bounds)
        @test !inBounds(T(max_x + eps(max_x), max_y), bounds)

        @test !inBounds(T(max_x, max_y + eps(max_y)), bounds)
    end

    # Test onbounds
    @testset "Testing onBounds()" for bounds in (Bounds{LLA{NAD27}}(0, 1, 78, 79),
                                                 Bounds{ENU}(-1, 1, -1, 1)
                                                )

        T = Geodesy.point_type(bounds)
        c = center(bounds)
        cx, cy = getX(c), getY(c)
        for _ = 1:1_000
            in_both =    T(cx + rand() - 0.5, cy + rand() - 0.5)
            in_x1 =      T(cx + rand() - 0.5, cy + rand() + 1.0)
            in_x2 =      T(cx + rand() - 0.5, cy - rand() - 1.0)
            in_y =       T(cx + rand() + 1.0, cy + rand() - 0.5)
            in_neither = T(cx + rand() + 1.0, cy + rand() + 1.0)

            @test onBounds(boundaryPoint(in_both, in_x1, bounds), bounds)
            @test onBounds(boundaryPoint(in_x2, in_both, bounds), bounds)
            @test onBounds(boundaryPoint(in_both, in_y, bounds), bounds)
            @test onBounds(boundaryPoint(in_neither, in_both, bounds), bounds)
        end
    end
end
