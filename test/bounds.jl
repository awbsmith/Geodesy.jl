using Geodesy
using Base.Test
using Geodesy: onBounds, inBounds, boundaryPoint, NAD27

# Construction

# inBounds

function test_bounds{T}(bounds::Bounds{T})

    min_x, min_y, max_x, max_y = bounds.min_x, bounds.min_y, bounds.max_x, bounds.max_y

    @test inBounds(T(min_x, min_y), bounds)
    @test !inBounds(T(min_x - eps(min_x), min_y), bounds)
    @test !inBounds(T(min_x, min_y - eps(min_y)), bounds)

    @test inBounds(T(max_x, max_y), bounds)
    @test !inBounds(T(max_x + eps(max_x), max_y), bounds)
    @test !inBounds(T(max_x, max_y + eps(max_y)), bounds)
end

for bounds in (Bounds{ENU}(1.1, 2.2, 3.3, 4.4),
               Bounds{LLA{NAD27}}(1.1, 2.2, 3.3, 4.4),
               )
    test_bounds(bounds)
end

# boundaryPoint

function test_boundary{T}(bounds::Bounds{T})
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

for bounds in (Bounds{LLA{NAD27}}(0, 1, 78, 79),
               Bounds{ENU}(-1, 1, -1, 1))
    test_boundary(bounds)
end
