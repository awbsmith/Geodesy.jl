using Base.Test
using Geodesy

################################################
### Helpers for testing approximate equality ###
################################################

# TODO: Move this to Compat.jl
if VERSION < v"0.4.0-dev+3616"
    fieldnames = names
end

macro tester(a, b)
    quote
        fieldnames($(esc(a)))[1]
    end
end

macro type_approx_eq(a, b)
    quote
        @test fieldnames($(esc(a))) == fieldnames($(esc(b)))
        for n in fieldnames($(esc(a)))
            @test getfield($(esc(a)), n) ≈ getfield($(esc(b)), n)
        end
    end
end

macro xyz_approx_eq(a, b)
    quote
        @test getX($(esc(a))) ≈ getX($(esc(b)))
        @test getY($(esc(a))) ≈ getY($(esc(b)))
        @test getZ($(esc(a))) ≈ getZ($(esc(b)))
    end
end
macro xy_approx_eq(a, b)
    quote
        @test getX($(esc(a))) ≈ getX($(esc(b)))
        @test getY($(esc(a))) ≈ getY($(esc(b)))
    end
end

macro xyz_approx_eq_eps(a, b, eps)
    quote
        @test ≈(getX($(esc(a))), getX($(esc(b))); atol=$(esc(eps)))
        @test ≈(getY($(esc(a))), getY($(esc(b))); atol=$(esc(eps)))
        @test ≈(getZ($(esc(a))), getZ($(esc(b))); atol=$(esc(eps)))
    end
end
macro xy_approx_eq_eps(a, b, eps)
    quote
        @test ≈(getX($(esc(a))), getX($(esc(b))); atol=$(esc(eps)))
        @test ≈(getY($(esc(a))), getY($(esc(b))); atol=$(esc(eps)))
    end
end
macro z_approx_eq_eps(a, b, eps)
    quote
        @test ≈(getZ($(esc(a))), getZ($(esc(b))); atol=$(esc(eps)))
    end
end

# GAH: Having issues passing the keyword "atol" into


# and run tests
for f in ["point", "bounds", "transform", "distance"]
    include("$f.jl")
end


#############################################
### Decimal <=> Degrees, Minutes, Seconds ###
#############################################

@testset "Testing DMS conversions" for (decimal, d, m, s) in [(0.013, 0.0, 0.0, 46.8),
                                                              (-0.013, -0.0, 0.0, 46.8),
                                                              (-0.263, -0.0, 15.0, 46.8),
                                                              (-179.51, -179.0, 30.0, 36.0)]
    @test Geodesy.dms2decimal(d, m, s) === decimal
    d2, m2, s2 = Geodesy.decimal2dms(decimal)
    @test d2 === d
    @test m2 === m
    @test s2 ≈ s
end
