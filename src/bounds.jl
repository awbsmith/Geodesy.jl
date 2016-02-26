
import Base.getindex, Base.setindex!, Base.*

###################
### Bounds Type ###
###################

type Bounds{T <: Geodesy_fam} 
	min_x::Float64
    max_x::Float64
    min_y::Float64
    max_y::Float64
end

# fill in the ellipse for the point type if its not provided
call(::Type{Bounds{LLA}}, min_x, max_x, min_y, max_y) = Bounds{add_param(LLA)}(min_x, max_x, min_y, max_y)
call(::Type{Bounds{LL}}, min_x, max_x, min_y, max_y) =  Bounds{add_param(LL)}(min_x, max_x, min_y, max_y)
call(::Type{Bounds{ECEF}}, min_x, max_x, min_y, max_y) = Bounds{add_param(ECEF)}(min_x, max_x, min_y, max_y)
call(::Type{Bounds{ENU}}, min_x, max_x, min_y, max_y) = Bounds{add_param(ENU)}(min_x, max_x, min_y, max_y)

# allow integer indexing
# TODO: learn the Julia convention for this, there must be a not horrible way...
Base.getindex{T}(bbox::Bounds{T}, idx::Int) = getfield(bbox, idx)
Base.setindex!{T}(bbox::Bounds{T}, val::Real, idx::Int) = setfield!(bbox, idx, Float64(val))  # why does Julia have reversed order for these?

# for degs / rads conversions
*{T}(bbox::Bounds{T}, val::Real) = Bounds{T}(bbox.min_x * val, bbox.max_x * val, bbox.min_y * val, bbox.max_y * val)


#############################
### Calculate LL bounds   ###
#############################

function call{T <: LL_fam}(::Type{Bounds{T}}, X::Vector{T}, degs::Bool=true)

	# build both paths around the circle to the first two points
	bounds_hyps = init_hyps(X, degs)

	# and convert both hypotheses
    pick_hyp(bounds_hyps, degs)

end

################################
### Calculate other bounds   ###
################################

function call{T <: Geodesy_fam}(::Type{Bounds{T}}, X::Vector{T}, degs::Bool=true)

	max_x = max_y = -Inf
    min_x = min_y = Inf
	for Xp in X
		min_x, max_x = min(getX(Xp), min_x), max(getX(Xp), max_x)
		min_y, max_y = min(getY(Xp), min_y), max(getY(Xp), max_y)
	end
	return Bounds{T}(min_x, max_x, min_y, max_y)
end

# unspecified type
call{T <: Geodesy_fam}(::Type{Bounds}, X::Vector{T}, degs::Bool=true) = Bounds{T}(X, degs)


################################
### Transform Bounds objects ###
################################

# there's not an unambiguous conversion, but for now,
# returning the minimum bounds that contain all points contained
# by the input bounds

# renamed it because it constructs a Bounds type not a ENU Type.  Does it belong as a transform?
function call{T <: Geodesy_fam, U <: LL_fam}(::Type{Bounds{T}}, bounds::Bounds{U}, ll_ref::U = center(bounds))

	oT = add_param(T)

    max_x = max_y = -Inf
    min_x = min_y = Inf

    xs = [bounds.min_x, bounds.max_x]
    ys = [bounds.min_y, bounds.max_y]
    
	# AWBS - what's this about?
	if bounds.min_y < 0.0 < bounds.max_y
        push!(ys, 0.0)
    end

    ref_x = getX(ll_ref)
    
	if bounds.min_x < ref_x < bounds.max_x ||
       (bounds.min_x > bounds.max_x && !(bounds.min_x >= ref_x >= bounds.max_x))
        push!(xs, ref_x)
    end

    for x_ll in xs, y_ll in ys
        pt = geotransform(oT, U(x_ll, y_ll), ll_ref)
        x, y = getX(pt), getY(pt)

        min_x, max_x = min(x, min_x), max(x, max_x)
        min_y, max_y = min(y, min_y), max(y, max_y)
    end

    return Bounds{oT}(min_x, max_x, min_y, max_y)
end


#############################
### Convenience Functions ###
#############################

### Get center point of a Bounds region ###
function center{T}(bounds::Bounds{T})
    x_mid = (bounds.min_x + bounds.max_x) / 2
    y_mid = (bounds.min_y + bounds.max_y) / 2
    return T(x_mid, y_mid)
end



### Check whether a location is within bounds ###
function inBounds{T <: Geodesy_fam}(loc::T, bounds::Bounds{T})
    x, y = getX(loc), getY(loc)

    bounds.min_x <= x <= bounds.max_x &&
    bounds.min_y <= y <= bounds.max_y
end


# only for points that have passed the inBounds test
function onBounds{T <: Union{LL, LLA, ENU}}(loc::T, bounds::Bounds{T})
    x, y = getX(loc), getY(loc)

    x == bounds.min_x || x == bounds.max_x ||
    y == bounds.min_y || y == bounds.max_y
end

# only for points where inBounds(p1) != inBounds(p2)
# TODO: return actual altitude rather than zero
function boundaryPoint{T <: Geodesy_fam}(p1::T, p2::T, bounds::Bounds{T})
    x1, y1 = getX(p1), getY(p1)
    x2, y2 = getX(p2), getY(p2)

    x, y = Inf, Inf

    # Move x to x bound if segment crosses boundary
    if x1 < bounds.min_x < x2 || x1 > bounds.min_x > x2
        x = bounds.min_x
        y = y1 + (y2 - y1) * (bounds.min_x - x1) / (x2 - x1)
    elseif x1 < bounds.max_x < x2 || x1 > bounds.max_x > x2
        x = bounds.max_x
        y = y1 + (y2 - y1) * (bounds.max_x - x1) / (x2 - x1)
    end

    p3 = T(x, y)
    inBounds(p3, bounds) && return p3

    # Move y to y bound if segment crosses boundary
    if y1 < bounds.min_y < y2 || y1 > bounds.min_y > y2
        x = x1 + (x2 - x1) * (bounds.min_y - y1) / (y2 - y1)
        y = bounds.min_y
    elseif y1 < bounds.max_y < y2 || y1 > bounds.max_y > y2
        x = x1 + (x2 - x1) * (bounds.max_y - y1) / (y2 - y1)
        y = bounds.max_y
    end

    p3 = T(x, y)
    inBounds(p3, bounds) && return p3

    error("Failed to find boundary point.")
end






#############################
### Helpers               ###
#############################

# bound theta in rads, so on output -pi <= theta < pi
bound_theta(theta::Real) = theta - floor((theta+pi) / (2*pi)) * 2*pi

# bound theta in degs, so on output -pi <= theta < pi
bound_thetad(theta::Real) = theta - floor((theta+180) / (360)) * 360

# function to start bounding boxes for the vector of LLA points X
function init_hyps{T <: Union{LLA, LL}}(X::Vector{T}, degs::Bool=true)

	bbox_hyps = [Bounds{T}(NaN, NaN, NaN, NaN), 
				 Bounds{T}(NaN, NaN, NaN, NaN)]  # two hypotheses
	sc = degs ? pi / 180 : 1.0

	# check for points
	if (length(X) > 0)
		# unsmart method - build both paths around the circle to the first two points
		cnt_min = length(X)+1
		for dim = 1:2

			# start point		
			ts = bound_theta(sc * (X[1][dim]))
			cnt = min(2,length(X))

			# find a second unique point
			while (cnt <= length(X)) && (ts == bound_theta(sc * X[cnt][dim]))
				cnt += 1
			end
			te = (cnt <= length(X)) ? bound_theta(sc * X[cnt][dim]) : ts + eps(ts)

			# build the two paths possible paths between points
			bbox_hyps[1][(dim-1)*2 + 1] = ts
			bbox_hyps[1][(dim-1)*2 + 2] = te > ts ? te : te + 2 * pi

			bbox_hyps[2][(dim-1)*2 + 1] = te
			bbox_hyps[2][(dim-1)*2 + 2] = te > ts ? ts + 2 * pi : ts

			cnt_min = min(cnt_min, cnt)  # where to start
		end

		if (degs)
			bbox_hyps[1] *= 180.0 / pi
			bbox_hyps[2] *= 180.0 / pi
		end

		# now include the rest of the points
		updatebounds!(bbox_hyps, X, cnt_min+1, degs)

	end

	return bbox_hyps # make a matrix out of it

end

# function to update a bounding box for an extra point
# updates bbox_h1 and bbox_h1
function updatebounds!{T <: Union{LLA, LL}}(bbox_hyps::Vector{Bounds{T}}, X::Vector{T}, first::Int=1, degs::Bool=true)

	sc = degs ? pi / 180 : 1.0
	bbox_hyps[1] *= sc
	bbox_hyps[2] *= sc

	# now add everything else top the competing paths
	for i = first:length(X)

		# x
		btheta = bound_theta(sc * getX(X[i]))
	    updatebounds_worker!(bbox_hyps[1], btheta, 1)
	    updatebounds_worker!(bbox_hyps[2], btheta, 1)

		# y
		btheta = bound_theta(sc * getY(X[i]))
	    updatebounds_worker!(bbox_hyps[1], btheta, 2)
	    updatebounds_worker!(bbox_hyps[2], btheta, 2)
	end

	# put it back in degs
	bbox_hyps[1] *= (1 / sc)
	bbox_hyps[2] *= (1 / sc)

end

# worker function for updating bounds
# input theta should be -pi <= theta < pi
function updatebounds_worker!{T <: Union{LLA, LL}}(bbox::Bounds{T}, theta::Real, dim::Int)

	# indexing shortcuts	
	idx1 = (dim-1)*2+1
	idx2 = (dim-1)*2+2

	# clockwise distance from the lower bound
    dlower = mod(theta - bbox[idx1], 2*pi)  # clockwise distance from lower bound

    # anticlockwise distance between bounds
    dbound = bbox[idx2] - bbox[idx1]

    if (dlower > dbound)

        # clockwise distance from the upper bound
        dupper = mod(theta - bbox[idx2], 2*pi)  

        # convert distance from the lower bound to anticlockwise
        dlower = 2*pi - dlower

        # pick the minimum adjustement to the bounds
        if (dlower < dupper)
            bbox[idx1] -= dlower
        else
            bbox[idx2] += dupper
        end
    end
end

# function to turn bounding boxes for each hypothesis into a Bounding box type
function pick_hyp{T <: Union{LLA, LL}}(bbox_hyps::Vector{Bounds{T}}, degs::Bool=true)

	wrap = degs ? 360.0 : 2*pi

	# x
	idx = 1 + (bbox_hyps[1].max_x - bbox_hyps[1].min_x >=  bbox_hyps[2].max_x - bbox_hyps[2].min_x)
	d1b = [bbox_hyps[idx].min_x, bbox_hyps[idx].max_x]
	d1b -= floor((mean(d1b) + wrap / 2) / (wrap)) * wrap  # keep it pretty

	# y
	idx = 1 + (bbox_hyps[1].max_y - bbox_hyps[1].min_y >=  bbox_hyps[2].max_y - bbox_hyps[2].min_y)
	d2b = [bbox_hyps[idx].min_y, bbox_hyps[idx].max_y]
	d2b -= floor((mean(d2b) + wrap / 2) / (wrap)) * wrap  # keep it pretty

	# and return the best combo
    return Bounds{T}(d1b[1], d1b[2], d2b[1], d2b[2])

end







