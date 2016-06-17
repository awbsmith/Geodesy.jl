import Base: convert, show

# datum traits and default values
datum_prop_fcns = [(:datum_dynamic, false),
                   (:datum_style, :(error("The datum style is unknown"))),
                   (:datum_ellipse, :(error("The ellipse is unknown"))),
                   (:datum_position, :(error("The reference position is unknown / undefined"))),
                   (:datum_time, :(error("The reference time is unknown / undefined")))
                  ]

# template cs properaties
function create_datum_templates()
    qb = quote; end
    for fcn in datum_prop_fcns
        push!(qb.args, :(@inline $(fcn[1]){T}(::Type{T}) = $(fcn[2])))
        push!(qb.args, :(@inline $(fcn[1]){T}(::T) = $(fcn)(T)))
    end
    return qb
end
create_datum_templates()


#
# define datums
#
abstract AbstractDatum
# show{T <: AbstractDatum}(io::Base.IO, ::Type{T}) = print(io,replace(string(T), "Geodesy2.", ""))

# return this for abstract datums
get_datum{T <: AbstractDatum}(::T) = UnknownDatum

#
# Unknown datum
#

immutable UnknownDatum <: AbstractDatum; end

datum_pos(::Type{UnknownDatum}) = error("Can't retrieve the reference position from an unknown datum")
datum_time(::Type{UnknownDatum}) = error("Can't retrieve the time of an unknown datum")
show(io::Base.IO, ::Type{UnknownDatum}) = print(io, "UnknownDatum")
show(io::Base.IO, ::UnknownDatum) = print(io, "UnknownDatum()")

# always allow conversion to an UnknownCS
convert(::Type{UnknownDatum}, datum::AbstractDatum) = UnknownDatum() # make this anything, or just AbstractDatums?



#
# Simple datums - the one datum defines all coordinates
#
abstract SimpleDatum <: AbstractDatum

abstract SimpleDynDatum <: SimpleDatum            # dynamic datums - need a time to transform between datums
@inline is_dynamic_datum{T <: SimpleDynDatum}(::Type{T}) = true


#
# Compund datums (e.g. one datum for horizontal coordinates and another one for vertical coordinates)
#
abstract CompoundDatum <: AbstractDatum
abstract CompoundDynDatum <: CompoundDatum   # dynamic datums - need a time to transform between datums

@inline is_dynamic_datum{T <: CompoundDynDatum}(::Type{T}) = true


#
# Horizontal / vertical compund datums are common so template for them
#
immutable HorzVertDatum{HDATUM <: SimpleDatum, VDATUM <: SimpleDatum} <: CompoundDatum; end
@inline is_dynamic_datum{H,V}(::Type{HorzVertDatum{H,V}}) = is_dynamic_datum(H)


#
# Create a template for geoid based vertical datums
#
abstract GeoidDatum <: SimpleDatum
@inline style{T <: GeoidDatum}(::Type{T}) = :GEOID


#
# Define a concrete datum type for when the datum is a point
#
immutable PositionDatum{POS} <: SimpleDatum;
    origin::POS
end
@inline style{T <: PositionDatum}(::Type{T}) = :REFERENCE_POSITION
@inline datum_position(datum::PositionDatum) = get_parameters(datum.origin)


#
# Define a concrete datum type for when the datum is a point
#
immutable TimedPositionDatum{TIME, POS} <: CompoundDatum;
    time::TIME
    origin::POS
end
@inline style{T <: TimedPositionDatum}(::Type{T}) = :TIMED_REFERENCE_POSITION
@inline datum_position(datum::TimedPositionDatum) = get_parameters(pdatum.origin)
@inline datum_time(datum::TimedPositionDatum) = get_parameters(pdatum.time)


#
# Overload for infering the datum when its not defined
#

@inline get_output_datums(::UnknownDatum, ::UnknownDatum) = error("No datum has been specified")
@inline get_output_datums(::UnknownDatum, datum) = (datum, datum)
@inline get_output_datums(datum, ::UnknownDatum) = (datum, datum)


# Completely generic
@inline get_output_datums(datum1, datum2) = (datum1, datum2)


#
# For assessing if datums match
#
@inline combine_datums(::UnknownDatum, ::UnknownDatum) = error("No datum has been specified")
@inline combine_datums(::UnknownDatum, datum) = datum
@inline combine_datums(datum, ::UnknownDatum) = datum

@inline combine_datums(::Type{UnknownDatum}, ::Type{UnknownDatum}) = UnknownDatum
@inline combine_datums{DATUM}(::Type{DATUM}, ::Type{UnknownDatum}) = DATUM
@inline combine_datums{DATUM}(::Type{UnknownDatum}, ::Type{DATUM}) = DATUM

@inline combine_datums{DATUM}(::Type{DATUM}, ::Type{DATUM}) = DATUM
@inline function combine_datums{DATUM}(d1::DATUM, d2::DATUM)
    (d1 == d2) || error("The Geodesy packaged doesn't perform datum transforms.  Please consider using Proj4 to handle datum transformations")
    return d1
end
@inline combine_datums(d1, d2) = error("The Geodesy packaged doesn't perform datum transforms.  Please consider using Proj4 to handle datum transformations")

