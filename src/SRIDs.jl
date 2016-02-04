#module SRIDs


# all epsg / esri comes from proj4
using Proj4
import Proj4.Projection  # going to overload this
import Proj4.transform  # going to overload this

abstract SRID_Types
abstract EPSG_Types <: SRID_Types # make epsg codes types
abstract ESRI_Types <: SRID_Types # make epsg codes types


# generate a symbol
macro srid_symgen(prefix, code)
	:(symbol($prefix*string($code)))
end


@doc """
	Create an SRID type from the given string, e.g.

	srid"EPSG4326"

""" ->
macro srid_str(text)

	# split into authority and code
	authority = UTF8String(matchall(r"\D+", text)[1])
	code = parse(matchall(r"\d+", text)[1])

	# create the new symbol
	authority = uppercase(authority)
	sym = @srid_symgen(authority, code)
	
	check_dict = authority == "EPSG" ? Proj4.epsg : authority == "ESRI" ? Proj4.esri :error("Unknown SRID Authority")
	if !haskey(check_dict, code)
		error("Proj4 does not support $(string(sym))")
	end

	# add the projection info
	proj_str = check_dict[code] 
	P4 = eval(:(Proj4.Projection($proj_str)))
	supertype = eval(:(symbol($authority*"_Types")))  

	# N.B. The quote block approach tries to add the new type to the Geodesey module
	#quote
	#	type $sym <: $supertype; end
	#	Proj4.Projection(::Type{$sym}) = $P4
	#end

	# The below expression statements add the new type to the current module instead
	exprs = Vector{Any}(0)
	expr = :(type $sym <: $supertype end)
	push!(exprs, expr)

	# make the type expression to get the Proj4 projection thingy
	expr = :(Proj4.Projection(::Type{$sym}) = $P4)
	push!(exprs, expr)
	return esc(Expr(:block, exprs...))

	
end


# define type for common use SRIDs now
srid"EPSG4326"    # initialise the wgs84 (GPS) EPSG type
typealias  EPSG_WGS84     EPSG4326

# build a Proj4 projection for each ellipse in out datum

"+proj=longlat +a=6378160 +b=6356774.50408554 +no_defs"


#############################################
# Calling the below macros to add all SRIDs 
# known to Proj4 kills Julia's typing system
#############################################

# macro to a type for each epsg code, as well as a function to get the code as an int 
macro populate_epsg()

	codes = keys(Proj4.epsg)
	codes = sort(collect(codes))

	# somewhere to put them all
	exprs = Vector{Any}()
	for code in codes
	
		# make the type expression
		epsg_sym = @srid_symgen("EPSG", eval(:($code))) #eval(:($code))????
		expr = :(type $epsg_sym <: EPSG_Types end)
		push!(exprs, expr)

		# make the type expression to get the Proj4 projection thingy
		str = Proj4.epsg[code] 
		expr = :(Proj4.Projection(::Type{$epsg_sym}) = Proj4.Projection($str))
		push!(exprs, expr)

	end
	return esc(Expr(:block, exprs...))
end


# macro to a type for each epsg code, as well as a function to get the code as an int 
macro populate_esri()

	codes = keys(Proj4.esri)
	codes = sort(collect(codes))

	# somewhere to put them all
	exprs = Vector{Any}()
	for code in codes
	
		# make the type expression
		esri_sym = @srid_symgen("ESRI", eval(:($code))) #eval(:($code))????
		expr = :(type $esri_sym <: ESRI_Types end)
		push!(exprs, expr)

		# make the type expression to get the Proj4 projection thingy
		str = Proj4.esri[code] 
		expr = :(Proj4.Projection(::Type{$esri_sym}) = Proj4.Projection($str))
		push!(exprs, expr)

	end
	return esc(Expr(:block, exprs...))
end


# and build
#@populate_epsg()
#@populate_esri()



#end

