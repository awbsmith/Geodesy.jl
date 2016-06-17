#
# Module to define tests related to CRS's
#
import Geodesy: UnknownCS, UnknownEllipse, eWGS84, WGS84_Ellipse

@testset "Datum tests" begin

    # check creating a dtum type
    @testset "Datum instantiation tests" begin

        # elliptic datums
        @test EllipticDatum(UnknownEllipse()) == EllipticDatum(UnknownEllipse)
        @test EllipticDatum(WGS84_Ellipse()) == EllipticDatum(WGS84_Ellipse)

        # position datums
        @test PositionDatum(UnknownPosition()) == PositionDatum(UnknownPosition)
    end

    @testset "Datum instantiation inference" begin

        # elliptic datums
        @inferred EllipticDatum(UnknownEllipse())
        @inferred EllipticDatum(UnknownEllipse)
        @inferred EllipticDatum(WGS84_Ellipse())
        @inferred EllipticDatum(WGS84_Ellipse)

        # position datums
        @inferred PositionDatum(UnknownPosition())
        @inferred PositionDatum(UnknownPosition)

        # unknown datums
        @inferred UnknownEllipticDatum()
        @inferred UnknownDatum()

        # value type datums
        @inferred EllipticDatum(eWGS84)

    end

    @testset "Unknown datum equality" begin

        # elliptic datums
        @test UnknownDatum() == EllipticDatum(UnknownEllipse())
        @test UnknownEllipticDatum() == EllipticDatum(UnknownEllipse())

        # position datums
        @test UnknownDatum() == PositionDatum(UnknownPosition())

    end

    @testset "Unknown datum equality inference" begin

        # elliptic datums
        @test UnknownEllipticDatum() == EllipticDatum(UnknownEllipse())

        # position datums
        @test UnknownDatum() == PositionDatum(UnknownPosition())

    end

    @testset "tagged vs value type equality" begin
        @test EllipticDatum(WGS84_Ellipse()) == EllipticDatum(eWGS84)
    end

end

@testset "CRS tests" begin

    # test unknowns
    @testset "CRS construction" begin

        @inferred CRS()
        @inferred CRS(UnknownCS(), UnknownDatum())
        @inferred UnknownCRS()

        @test CRS() == CRS(UnknownCS(), UnknownDatum())
        @test UnknownCRS() == CRS(UnknownCS(), UnknownDatum())
    end

    @testset "CRS datum promotion" begin

        # test promoting ellipses to elliptic datums
        @inferred CRS(UnknownCS(), EllipticDatum(eWGS84))
        @inferred CRS(UnknownCS(), eWGS84)
        @test CRS(UnknownCS(), EllipticDatum(eWGS84)) == CRS(UnknownCS(), eWGS84)

        @inferred CRS(UnknownCS(), EllipticDatum(WGS84_Ellipse))
        @inferred CRS(UnknownCS(), WGS84_Ellipse)
        @test CRS(UnknownCS(), EllipticDatum(WGS84_Ellipse)) == CRS(UnknownCS(), WGS84_Ellipse)

        # test promoting positions to position datums
        @inferred CRS(UnknownCS(), PositionDatum(GeoPosition((0.0,0.0,0.0))))
        @inferred CRS(UnknownCS(), GeoPosition((0.0,0.0,0.0)))
        @test CRS(UnknownCS(), PositionDatum(GeoPosition((0.0,0.0,0.0)))) == CRS(UnknownCS(), GeoPosition((0.0,0.0,0.0)))

        @inferred CRS(UnknownCS(), PositionDatum(LLA(0.0,0.0,0.0)))
        @inferred CRS(UnknownCS(), LLA(0.0,0.0,0.0))
        @test CRS(UnknownCS(), PositionDatum(LLA(0.0,0.0,0.0))) == CRS(UnknownCS(), LLA(0.0,0.0,0.0))
    end

end

