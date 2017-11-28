using Geodesy
using Base.Test

@testset "Testing Points" begin

    @testset "Testing Construction" begin

        x, y = (rand(2) - .5) * 10_000

        @test LLA(x, y) == LLA(x, y, 0.0)
        @test LLA(x, y) == LLA(x, y, 0)
        @test LLA_WGS84(x, y) == LLA_WGS84(x, y, 0)

        @test LL(x, y, 0.) == LL(x, y)
        @test LL(x, y, 0) == LL(x, y)
        @test LL{WGS84}(x, y, 0) == LL{WGS84}(x, y)

        lla = LLA(x, y, rand())
        @test ENU(x, y) == ENU(x, y, 0)
        @test ENU(x, y) == ENU(x, y, 0.0)
        @test ENU{lla}(x, y) == ENU{lla}(x, y, 0.0)


        ll = LL(x, y)
        @test LLA(getX(lla), getY(lla), getZ(lla)) == lla
        @test getY(ll) == y
        @test getX(ll) == x

        enu = ENU(x, y, rand())
        @test ENU(getX(enu), getY(enu), getZ(enu)) == enu

    end


    ################################################
    ### Testing parameter stripping and addition ###
    ################################################

    @testset "Testing Parameter Adding / Removal" begin

        lla_wgs84 = LLA_WGS84(-71.0960, 42.3673, 0)
        lla_null = LLA_NULL(-71.0960, 42.3673, 0)
        @xyz_approx_eq convert(LLA_NULL, lla_wgs84) lla_null
        @xyz_approx_eq convert(LLA{WGS84}, lla_null) lla_wgs84

    end

end
