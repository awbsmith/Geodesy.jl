
################################################
### Helpers for testing approximate equality ###
################################################

macro approx_eq(a, b)
    quote
        @test Geodesy2.numel($(a)) == Geodesy2.numel($(b))
        for i = 1:Geodesy2.numel($(a))
            @test_approx_eq $(a)[i] $(b)[i]
        end
    end
end

macro approx_eq_eps(a, b, eps)
    quote
        @test Geodesy2.numel($(a)) == Geodesy2.numel($(b))
        for i = 1:Geodesy2.numel($(a))
            @test abs($(a)[i] - $(b)[i]) < $(eps)
        end
    end
end


###################################
### Testing fixed relationships ###
###################################



@testset "Point Transforms" begin

    # these should be a pair
    test_points = [Dict(:NAME => "(no height)",
                        :LLA => LLA{WGS84}(-71.0960, 42.3673, 0),
                        :ECEF => ECEF{WGS84}(1529073.1560519305057823658, -4465040.0190131030976772308, 4275835.3392603071406483650),  # should match LLA
                        :ENU => ENU(-343.4937490839769793638, 478.7648554667879920999, -0.0272428852242967423),                        # matches LLA when centered on LLA_REF
                        :LLA_REF => LLA{WGS84}(-71.09183, 42.36299, 0.0)

                       ),

                   # add height to this one
                   Dict(:NAME => "(with height)",
                        :LLA => LLA{WGS84}(-71.0960, 42.3673, 10.0),
                        :ECEF => ECEF{WGS84}(1529075.5497715696692466736, -4465047.0089036403223872185, 4275842.0780685534700751305),  # should match LLA
                        :ENU => ENU(-343.4942868136923834754, 478.7656077168462047666, 4.9727570727289105434),                         # matches LLA when centered on LLA_REF
                        :LLA_REF => LLA{WGS84}(-71.09183, 42.36299, 5.0)
                       )

                  ]


    @testset "Transform type: $(test_point[:NAME])" for test_point in test_points

        # LLA -> ECEF
        @testset "LLA -> ECEF" begin
            ecef_test = ECEF(test_point[:LLA])
            @approx_eq_eps ecef_test test_point[:ECEF] 1e-8
            @test typeof(ecef_test) == typeof(test_point[:ECEF])  # check the datum inference as well
        end

        # ECEF -> LLA
        @testset "ECEF -> LLA" begin
            lla_test = LLA(test_point[:ECEF])
            @approx_eq_eps lla_test test_point[:LLA]  1e-8
            @test typeof(lla_test) == typeof(test_point[:LLA])  # check the datum inference as well
        end

        # LLA -> ENU
        @testset "LLA -> ENU" begin
            enu_test = ENU(test_point[:LLA], test_point[:LLA_REF])
            @approx_eq_eps enu_test test_point[:ENU] 1e-8
            @test typeof(enu_test) == typeof(test_point[:ENU])   # check the datum inference as well
        end

        # ECEF -> ENU
        @testset "ECEF -> ENU" begin
            enu_test = ENU(test_point[:ECEF], test_point[:LLA_REF])
            @approx_eq_eps enu_test test_point[:ENU] 1e-8
            @test typeof(enu_test) == typeof(test_point[:ENU])  # check the datum inference as well
        end

        # ENU -> LLA
        @testset "ENU -> LLA" begin
            lla_test = LLA(test_point[:ENU], test_point[:LLA_REF])
            @approx_eq_eps lla_test  test_point[:LLA] 1e-8
            @test typeof(lla_test) == typeof(test_point[:LLA])  # check the datum inference as well
        end

        # ENU -> ECEF
        @testset "ENU -> ECEF" begin
            ecef_test = ECEF(test_point[:ENU], test_point[:LLA_REF])
            @approx_eq_eps ecef_test  test_point[:ECEF] 1e-8
            @test typeof(ecef_test) == typeof(test_point[:ECEF])  # check the datum inference as well
        end
    end
end

