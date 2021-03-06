using DataFrames
using StatsBase


"""
	in_range(x, xmin, xmax)

Given vector x, select values between xmin and xmax
"""
function in_range(x::Vector{T}, xmin::T, xmax::T) where T <: Number
    xi = x[x.>xmin]
    xu = xi[xi.<xmax]
    return xu
end


"""
	rxy(x::T,y::T) where T
	r -> sqrt{x^2 + y^2}
"""
function rxy(x::T,y::T) where T <: Number
    return sqrt(x^2 + y^2)
end


"""
	hixy(x::T,y::T) where T
	phi -> atan(y/x)
"""
function phixy(x::T,y::T) where T <: Number
    return atan(y,x)
end

"""
	phixy(hitdf::DataFrame)
	phi -> atan(y/x) where y and x are columns of the data frame
"""
function fphi(hitdf::DataFrame)
    return atan.(hitdf.y,hitdf.x)
end


"""
	dxyz(x1::Vector{T}, x2::Vector{T})

Distance between two points.
"""

function dxyz(x1::Vector{T}, x2::Vector{T}) where T <: Number
    return sqrt((x1[1] - x2[1])^2 + (x1[2] - x2[2])^2 + (x1[3] - x2[3])^2)
end

"""
function wstd(x::Vector{T}, q::Vector{T}) where T

Compute the std deviation in x weighted by q:
Sqrt(1/Q Sum_i (x - x_mean) * qi )
"""
function wstd(x::Vector{T}, q::Vector{T}) where T <: Number
	xmean = mean(x)
	qs = sum((x.-xmean).^2 .* q)
	Q = sum(q)
	return sqrt(qs/Q)
end
#wstd(x::Vector{Float32}, q::Vector{Float32}) = gwstd(x, q)
#wstd(x::Vector{Float64}, q::Vector{Float64}) = gwstd(x, q)
#wstd(x::Vector{Number}, q::Vector{Number}) = gwstd(x, q)


"""
	mean_std(x, xmin, xmax)
	Returns mean and std for a vector x in the interval between xmin and xmax
"""
function mean_std(x::Vector{T}, xmin::T, xmax::T) where T <: Number
    xx = in_range(x, xmin, xmax)
    xm = mean(xx)
    xs = std(xx)
    return xm, xs
end
#mean_std(x::Vector{Float32}, xmin::Float32, xmax::Float32) = gmean_std(x, xmin, xmax)
#mean_std(x::Vector{Float64}, xmin::Float64, xmax::Float64) = gmean_std(x, xmin, xmax)
#mean_std(x::Vector{Int64}, xmin::Int64, xmax::Int64) = gmean_std(x, xmin, xmax)
