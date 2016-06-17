# for different ways to instantiate stuffs
@inline get_parameters(data) = data                                  # usual
@inline get_parameters{DTYPE}(::Type{DTYPE}) = DTYPE()               # overload this!
@inline get_parameters{POS}(::Type{Val{POS}}) = get_parameters(POS)  # strip the Val and go again
@inline get_parameters{POS}(::Val{POS}) = get_parameters(POS)        # the data is in POS


#############################################
### Decimal <=> Degrees, Minutes, Seconds ###
#############################################

function decimal2dms(x::Float64)
    d = trunc(x, 0)
    ms = 60 * abs(x - d)
    m = trunc(ms, 0)
    s = 60 * rem(ms, 1)

    return d, m, s
end

function dms2decimal(d::Float64, m::Float64, s::Float64)
    signbit(d) ? d - m/60 - s/3600 :
                 d + m/60 + s/3600
end
