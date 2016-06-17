
#############################################
### Decimal <=> Degrees, Minutes, Seconds ###
#############################################

@testset "Utility Tests" begin
	@testset "Lat / Lon format conversions" for (decimal, d, m, s) in [(0.013, 0.0, 0.0, 46.8),
				                   									   (-0.013, -0.0, 0.0, 46.8),
													                   (-0.263, -0.0, 15.0, 46.8),
													                   (-179.51, -179.0, 30.0, 36.0)]

		@test Geodesy.dms2decimal(d, m, s) === decimal
		d2, m2, s2 = Geodesy.decimal2dms(decimal)
		@test d2 === d
		@test m2 === m
		@test_approx_eq s2 s
	end
end

