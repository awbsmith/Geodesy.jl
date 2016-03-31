import Base: (+), (-), isnan, (*)

#####################################################
# add basic maths to the point type
#####################################################

function add_maths{geodesy_type}(::Type{geodesy_type}, fields)

    nf = length(fields)  # number of fields and their names for the Geodesy variables

    # allow interaction with these types
    interact_types = [Vec{nf, Float64}, Vector]  

    # create an expression to add to
    qb = quote end

    #
    # Allow subtracting the same Geodesy representation
    #  
    rhs_expr = :(())
    append!(rhs_expr.args, [:(X.$(field) - Y.$(field)) for field in fields[1:nf]])

    # form the whole expression
    qn = quote
        (-){T <: $(geodesy_type)}(X::T, Y::T) = Vec{$(nf), Float64}($(rhs_expr.args...))
    end
    append!(qb.args, qn.args)

    #
    # Allow Rotations
    #
    if (nf == 3)
        qn = quote
            (*)(R::Mat{3,3,Float64}, X::$(geodesy_type)) = R * Vec(X)
        end
        append!(qb.args, qn.args)  # include them
    end

    #
    # Allow uniform scaling
    #
    # build the rhs epxression
    rhs_expr = :(())
    append!(rhs_expr.args, [:(s * X.$(field)) for field in fields[1:nf]])

    # form the whole expression
    qn = quote
        (*){T <: $(geodesy_type), U <: Real}(X::T, s::U) = T($(rhs_expr.args...))
        (*){T <: $(geodesy_type), U <: Real}(s::U, X::T) = T($(rhs_expr.args...))
    end
    append!(qb.args, qn.args)


    # nan check
    rhs_nan_expr = Expr(:||, :(isnan(X.$(fields[nf - 1]))), :(isnan(X.$(fields[nf]))))
    for i in nf-2:-1:1
        rhs_nan_expr = Expr(:||,  :(isnan(X.$(fields[i]))) , rhs_nan_expr)
    end
    

    #
    #  Fill function for the interactions types
    #
    for iT in interact_types

        # add a vector to it        
        rhs_add_expr = :(())
        append!(rhs_add_expr.args, [:(X.$(field) + dX[$(i)]) for (i, field) in enumerate(fields)])

        # subtract a vector from it        
        rhs_sub_expr = :(())
        append!(rhs_sub_expr.args, [:(X.$(field) - dX[$(i)]) for (i, field) in enumerate(fields)])

        qn = quote
            
            # addition 
            (+){T <: $(geodesy_type)}(X::T, dX::$(iT)) = T($(rhs_add_expr.args...))
            (+){T <: $(geodesy_type)}(dX::$(iT), X::T) = T($(rhs_add_expr.args...))

            # subtraction
            (-){T <: $(geodesy_type)}(X::T, dX::$(iT)) = T($(rhs_sub_expr.args...))

            # NaN
            isnan(X::$(geodesy_type)) = $(rhs_nan_expr)
            
        
        end
        append!(qb.args, qn.args)  # include them
    end
   
    
    return qb

end
