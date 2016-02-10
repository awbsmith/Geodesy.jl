
# auth is a symbol, code is an integer
immutable SRID{auth, code} <: Datum  end  # good to have it as a subtype of datum? 
show{auth, code}(io::IO, ::Type{SRID{auth, code}}) = print(io, "$(auth)$(code)")

# calling this something different to not overload Proj4 stuff with generated functions
@generated function get_projection{auth, code}(::Type{SRID{auth, code}})

	println("Gen: $(SRID{auth, code})")
	dict_sym = symbol(lowercase(string(auth)))

	
	local dict
	try # hasfield / isfield? 
		dict = Proj4.(dict_sym)
	catch
		error("Proj4 does not know the SRID Authority: $(auth)")
	end

	if !haskey(dict, code)
		error("Proj4 does not know the code $(code) for authority $(auth)")
	end

	# add the projection info
	proj = Proj4.Projection(dict[code])
	return :($proj)

end



