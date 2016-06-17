using BaseTestNext
using FixedSizeArrays
using Geodesy

@testset "Geodesy Tests" begin
    for f in ["crs", "position"] #, "misc", "transform"]
        include("$f.jl")
    end
end
