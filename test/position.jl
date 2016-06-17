import Geodesy2: ndims, default_coord_system
import Geodesy2: LLA_CS, LL_CS, ECEF_CS, ENU_CS, UnknownCS, AbstractPosition
import Geodesy2: eWGS84



@testset "Position Type Tests" begin

    #
    #  Test fixed relationships to check they aren't broken
    #
    @testset "Default Coordinate Systems" begin

        # test the CS is correct for each type
        @test default_coord_system(GeoPosition) == UnknownCS
        @test default_coord_system(LLA) == LLA_CS
        @test default_coord_system(LL) == LL_CS
        @test default_coord_system(ECEF) == ECEF_CS
        @test default_coord_system(ENU) == ENU_CS
    end

    #
    #  Now test constructors
    #

    @testset "$(pType)" for pType in setdiff(subtypes(AbstractPosition), [UnknownPosition])

        # the expected CRS for this type
        default_crs = CRS{default_coord_system(pType), UnknownDatum}

        # build one with no CRS info
        xv = randn(ndims(pType))
        X = pType(xv..., default_crs()) # build it using element by element construction

        @testset "Constructors" begin

            # test the crs
            @test get_crs(X) == default_crs()
            @test get_coord_system(X) == default_coord_system(pType)()
            @test get_datum(X) == UnknownDatum()

            # create a tuple and a FixedVector form

            xt = (xv...)
            xfv = Vec(xv)

            # try the various construction methods
            @testset "$(typeof(data))" for data in (xt, xv, xfv)

                # check CRS auto fill constructor
                @test X == pType(data)

                # check CRS instantiation
                @test X == pType(data, default_crs)

                # check CRS in the type
                @test X == pType{default_crs}(data)

                # check CRS and other parameters in the type
                @test X == pType{default_crs}(data)

                # this version is never going to be type safe when a type parameter comes from a variable length vector (manually specfify the size)
                if (pType != GeoPosition) || !isa(data, AbstractVector)
                    @inferred pType(data)
                    @inferred pType(data, default_crs)
                    @inferred pType{default_crs}(data)
                    @inferred pType{default_crs}(data)
                end
            end
        end
    end
end

