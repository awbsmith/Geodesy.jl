
###################
### Bounds Type ###
###################

type Bounds{T <: Union{LL, LLA, ENU}}
    min_y::Float64
    max_y::Float64
    min_x::Float64
    max_x::Float64
end


function Bounds(min_lat, max_lat, min_lon, max_lon)
	
	#= Removed -> I like to express min < max, and the lon range 179 -> 181 will fail this check
    if !(-90 <= min_lat <= max_lat <= 90 &&
         -180 <= min_lon <= 180 &&
         -180 <= max_lon <= 180)
        throw(ArgumentError("Bounds out of range of LL coordinate system. " *
                            "Perhaps you're looking for Bounds{ENU}(...)"))
    end
	=#
    Bounds{LL{WGS84_ELLIPSE}}(min_lat, max_lat, min_lon, max_lon)
end


#############################
### Calculate LL bounds   ###
#############################

" Functions to calculate the bounds of a vector of points "
function bounds{T <: Union{LLA, LL}}(X::Vector{T}, degs::Bool=true)

	# unsmart method - build both paths around the circle to the first two points
	(bbox_h1, bbox_h2) = init_bounds(X, degs)

	# and convert both hypotheses
	bounds_from_bboxes(bbox_h1, bbox_h2, degs)

end


#############################
### Convenience Functions ###
#############################

### Get center point of Bounds region ###
function center(bounds::Bounds{ENU})
    x_mid = (bounds.min_x + bounds.max_x) / 2
    y_mid = (bounds.min_y + bounds.max_y) / 2

    return ENU(x_mid, y_mid)
end

function center{T <: Union{LL, LLA}}(bounds::Bounds{T})
    x_mid = (bounds.min_x + bounds.max_x) / 2
    y_mid = (bounds.min_y + bounds.max_y) / 2

    if bounds.min_x > bounds.max_x
        x_mid = x_mid > 0 ? x_mid - 180 : x_mid + 180
    end

    return T(y_mid, x_mid)
end

### Check whether a location is within bounds ###
function inBounds(loc::ENU, bounds::Bounds{ENU})
    x, y = getX(loc), getY(loc)

    bounds.min_x <= x <= bounds.max_x &&
    bounds.min_y <= y <= bounds.max_y
end

function inBounds{T <: Union{LL, LLA}}(loc::T, bounds::Bounds{T})
    x, y = getX(loc), getY(loc)

    min_x, max_x = bounds.min_x, bounds.max_x

    (min_x > max_x ? !(max_x < x < min_x) : min_x <= x <= max_x) &&
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
function boundaryPoint{T <: Union{LL, LLA, ENU}}(p1::T, p2::T, bounds::Bounds{T})
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

    p3 = T(XY(x, y))
    inBounds(p3, bounds) && return p3

    # Move y to y bound if segment crosses boundary
    if y1 < bounds.min_y < y2 || y1 > bounds.min_y > y2
        x = x1 + (x2 - x1) * (bounds.min_y - y1) / (y2 - y1)
        y = bounds.min_y
    elseif y1 < bounds.max_y < y2 || y1 > bounds.max_y > y2
        x = x1 + (x2 - x1) * (bounds.max_y - y1) / (y2 - y1)
        y = bounds.max_y
    end

    p3 = T(XY(x, y))
    inBounds(p3, bounds) && return p3

    error("Failed to find boundary point.")
end

#############################
### Helpers               ###
#############################

# bound theta in degs, so on output -pi <= theta < pi
bound_theta(theta::Real) = theta - floor((theta+pi) / (2*pi)) * 2*pi

# function to start bounding boxes for the vector of LLA points X
function init_bounds{T <: Union{LLA, LL}}(X::Vector{T}, degs::Bool=true)

	bbox_h1 = fill(NaN, 4) #xmin, xmax, ymin, ymax
	bbox_h2 = fill(NaN, 4) #xmin, xmax, ymin, ymax
	sc = degs ? pi / 180 : 1.0

	# check for points
	if (length(X) > 0)

		# unsmart method - build both paths around the circle to the first two points
		cnt_min = length(X)+1
		for dim = 1:2

			# start point		
			ts = bound_theta(sc * X[1][dim])
			cnt = min(2,length(X))

			# find a second unique point
			while (cnt <= length(X)) && (ts == bound_theta(sc * X[cnt][dim]))
				cnt += 1
			end
			te = (cnt <= length(X)) ? bound_theta(sc * X[cnt][dim]) : ts + eps(ts)

			# build the two paths possible paths between points
			bbox_h1[(dim-1)*2+1] = ts
			bbox_h1[(dim-1)*2+2] = te > ts ? te : te+2*pi

			bbox_h2[(dim-1)*2+1] = te
			bbox_h2[(dim-1)*2+2] = te > ts ? ts+2*pi : ts

			cnt_min = min(cnt_min, cnt)  # where to start
		end

		# now include the rest of the points
		updatebounds!(bbox_h1, bbox_h2, X, cnt_min+1, degs)

	end


	return (bbox_h1, bbox_h2)

end

# function to update a bounding box for an extra point
# updates bbox_h1 and bbox_h1
function updatebounds!{T <: Union{LLA, LL}}(bbox_h1::AbstractVector, bbox_h2::AbstractVector, X::Vector{T}, first::Int=1, degs::Bool=true)

	sc = degs ? pi / 180 : 1.0

	# now add everything else top the competing paths
	for i = first:length(X)

		# first coord
		btheta = bound_theta(sc * X[i][1])
        updatebounds_worker!(btheta, 1, bbox_h1)
        updatebounds_worker!(btheta, 1, bbox_h2)

		# second coord
		btheta = bound_theta(sc * X[i][2])
        updatebounds_worker!(btheta, 2, bbox_h1)
        updatebounds_worker!(btheta, 2, bbox_h2)

	end
end

# worker function for updating bounds
# input theta should be -pi <= theta < pi
function updatebounds_worker!(theta::Real, dim::Int, bbox::AbstractVector)

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
function  bounds_from_bboxes(bbox_h1::AbstractVector, bbox_h2::AbstractVector, degs::Bool=true)

	sc = degs ? 180/pi : 1

	# first dim
	d1b = ((bbox_h1[2] - bbox_h1[1]) <= (bbox_h2[2] - bbox_h2[1]) ? bbox_h1[1:2] : bbox_h2[1:2]) 
	d1b -= floor((mean(d1b)+pi) / (2*pi)) * 2*pi  # keep it pretty
	d1b *= sc

	# second dim
	d2b = ((bbox_h1[4] - bbox_h1[3]) <= (bbox_h2[4] - bbox_h2[3]) ? bbox_h1[3:4] : bbox_h2[3:4]) 
	d2b -= floor((mean(d2b)+pi) / (2*pi)) * 2*pi  # keep it pretty
	d2b *= sc

	# and store the best one
    return Bounds(d1b[1], d1b[2], d2b[1], d2b[2])

end





