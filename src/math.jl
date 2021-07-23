using LinearAlgebra
using GLM
using LsqFit
using Distributions
using StatsBase

#linear fit wrapper

function lfit(ndf::DataFrame)
    lr = lm(@formula(y_mean ~ x_mean), ndf)
    c = coef(lr)
    return x -> c[1] + c[2]*x, predict(lr), c
end

function fit_pol2(x,y)
    @. pol(x, p) = p[1] + p[2] * x + p[3] * x^2
    p0 = [1.0, 1.0, 1.0]
    fq = curve_fit(pol, x, y, p0)
    cfq = coef(fq)
    @info "coef(fq)" cfq
    sfq = stderror(fq)
    @info "std(fq)" sfq
    @info "margin_of_error (90%)" margin_error(fq, 0.1)
    @info " confidence_interval (90%)" confidence_interval(fq, 0.1)
    return cfq
end


# function gfit_gauss(x, xmin, xmax, bins=50)
#     function gausx(x, μ, σ, N)
#         return N * exp(-(x - μ)^2/(2*σ^2))
#     end
#
#     xmu, xstd = mean_std(x, xmin, xmax)
#     @debug xmu xstd
#
#     h = hist1d(x,  bins, xmin, xmax)
#     edges = collect(h.edges[1])
#     w = h.weights
#     c =[0.5 *(edges[i] + edges[i+1]) for i in 1:length(edges)-1]
#     @debug "histo" edges w c
#
#     @. gauss1(x, p) = p[1]* exp(-(x - xmu)^2/(2*xstd^2))
#     p0 = [1.0]
#     fq = curve_fit(gauss1, c, w, p0)
#     NN =coef(fq)[1][1]
#     @debug "gauss1" NN
#
#     @. gauss3(x, p) = p[1]* exp(-(x - p[2])^2/(2*p[3]^2))
#     p0 = [NN, xmu, xstd]
#     fq = curve_fit(gauss3, c, w, p0)
#     cfq = coef(fq)
#     @debug "coef(fq)" cfq
#     NN = cfq[1]
#     mu =cfq[2]
#     std  =cfq[3]
#     @debug "gauss3" NN mu std
#
#     #return xmu, xstd, mu, std, NN, c, gausx.(c, (mu,), (std,), (NN),)
# 	return FGauss(mu, std, NN, h, c, gausx.(c, (mu,), (std,), (NN),))
# end
#fit_gauss(x::Vector{Float64},
#          xmin::Float64, xmax::Float64, bins::Integer=50) = gfit_gauss(x, xmin, xmax, bins)
#fit_gauss(x::Vector{Float32},
#          xmin::Float32, xmax::Float32, bins::Integer=50) = gfit_gauss(x, xmin, xmax, bins)

"""
	fit_gauss(y, xmin, xmax, bins=25)

Fit a normal distribution to data
"""

struct FGauss
	mu::Vector{Number}
	std::Vector{Number}
	C::Vector{Number}
	h::Histogram
	X::Vector{Number}
	Y::Vector{Number}
	g::Vector{Function}
end

function gausg(μ, σ, C)
	function gausx(x)
		return C * pdf(Normal(μ, σ,), x)
	end
	return gausx
end

function gausg2(μ1, σ1, C1, μ2, σ2, C2)
	function gausx(x)
		return C1 * pdf(Normal(μ1, σ1,), x) + C2 * pdf(Normal(μ2, σ2,), x)
	end
	return gausx
end

@. gauss1(x, p) = p[1]* pdf(Normal(p[2], p[3]), x)
@. gauss2(x, p) = p[1]* pdf(Normal(p[2], p[3]), x) +  p[4]* pdf(Normal(p[5], p[6]), x)


function cfit(ffit, x, y, p0, lb, ub)
	fq = curve_fit(ffit, x, y, p0, lower=lb, upper=ub)
    cfq = coef(fq)
    @debug "coef(fq)" cfq
    return cfq
end

function hfit_gauss(h::Histogram)
	c = centers(h)
	w = h.weights
	@debug "histo"  w c
	mu, sigma = mean_and_std(c, Weights(h.weights); corrected = false)
	@debug "mu, std" mu, sigma

	# fit parameters lb, ub, po are lower, upper bounds and pars
    lb = [1., mu - 5*sigma, 0.1]
    ub = [sum(w), mu + 5*sigma, 5*sigma]
    p0_bounds = [100., mu, sigma]
	CC, μ, σ  = cfit(gauss1, c, w, p0_bounds, lb, ub)
	gx = gausg(μ, σ, CC)
	return FGauss([μ], [σ], [CC], h, c, gx.(c), [gx])

end

function gfit_gauss(y, xmin, xmax, bins=25)

	# fit the unbinned distribution
	x  =  in_range(y, xmin, xmax)
	ft =  fit(Normal,x)
	@debug ft.μ ft.σ

	# bin distribution
    h = hist1d(x, bins, xmin, xmax)
    c = centers(h)
    w = h.weights
    @debug "histo"  w c

	# fit parameters lb, ub, po are lower, upper bounds and pars
    lb = [1., ft.μ - 5*ft.σ, 0.1]
    ub = [sum(w), ft.μ + 5*ft.σ, 5*ft.σ]
    p0_bounds = [100., ft.μ, ft.σ]
	CC, mu, sigma  = cfit(gauss1, c, w, p0_bounds, lb, ub)
	gx = gausg(mu, sigma, CC)

	return FGauss([mu], [sigma], [CC], h, c, gx.(c), [gx])
end


function gfit_gauss2(y, xmin, xmax, bins=25)

	# fit the unbinned distribution
	x  =  in_range(y, xmin, xmax)
	mu, sigma =  mean_and_std(x)
	@debug " mu sigma " mu sigma

	# bin distribution
    h = hist1d(x, bins, xmin, xmax)
    c = centers(h)
    w = h.weights
    @debug "histo"  w c

	# fit parameters lb, ub, po are lower, upper bounds and pars

    p0 = [sum(w)/10.,  mu,           sigma,   sum(w)/10.,  mu,           sigma]
	lb = [sum(w)/100., mu - 5*sigma, 0.1,     sum(w)/100., mu - 5*sigma, 0.1]
    ub = [sum(w),      mu + 5*sigma, 5*sigma, sum(w),      mu + 5*sigma, 5*sigma]

	fq = curve_fit(gauss2, c, w, p0, lower=lb, upper=ub)
    CC1, mu1, sigma1, CC2, mu2, sigma2   = coef(fq)
    @debug "coef(fq)" cfq
	gx = gausg2(mu1, sigma1, CC1, mu2, sigma2, CC2)
	gx1 = gausg(mu1, sigma1, CC1)
	gx2 = gausg(mu2, sigma2, CC2)

	return FGauss([mu1, mu2], [sigma1, sigma2], [CC1, CC2], h, c, gx.(c), [gx1, gx2])
end

fit_gauss(x::Vector{Float64},
        xmin::Float64, xmax::Float64, bins::Integer=50) = gfit_gauss(x, xmin, xmax, bins)
fit_gauss(x::Vector{Float32},
        xmin::Float32, xmax::Float32, bins::Integer=50) = gfit_gauss(x, xmin, xmax, bins)
fit_gauss2(x::Vector{Float64},
        xmin::Float64, xmax::Float64, bins::Integer=50) = gfit_gauss2(x, xmin, xmax, bins)
fit_gauss2(x::Vector{Float32},
        xmin::Float32, xmax::Float32, bins::Integer=50) = gfit_gauss2(x, xmin, xmax, bins)
fit_gauss(h::Histogram) = hfit_gauss(h::Histogram)
#misc math

"""
	rxy(x::Number,y::Number)
	r -> sqrt{x^2 + y^2}
"""
function rxy(x::Number,y::Number)
    return sqrt(x^2 + y^2)
end


"""
	phixy(x::Number,y::Number)
	phi -> atan(y/x)
"""
function phixy(x::Number,y::Number)
    return atan(y,x)
end

"""
	phixy(x::Number,y::Number)
	phi -> atan(y/x)
"""
function fphi(hitdf::DataFrame)
    return atan.(hitdf.y,hitdf.x)
end


"""
	dxyz(x1::Vector{Number}, x2::Vector{Number})

Distance between two points.
"""
function gdxyz(x1, x2)
    return sqrt((x1[1] - x2[1])^2 + (x1[2] - x2[2])^2 + (x1[3] - x2[3])^2)
end


function dxyz(x1::Vector{Float32}, x2::Vector{Float32})
    return gdxyz(x1, x2)
end
function dxyz(x1::Vector{Float64}, x2::Vector{Float64})
    return gdxyz(x1, x2)
end
function dxyz(x1::Vector{Int64}, x2::Vector{Int64})
    return gdxyz(x1, x2)
end


"""
function gwstd(x, q)

Compute the std deviation in x weighted by q:
Sqrt(1/Q Sum_i (x - x_mean) * qi )
"""
function gwstd(x, q)
	xmean = mean(x)
	qs = sum((x.-xmean).^2 .* q)
	Q = sum(q)
	return sqrt(qs/Q)
end
wstd(x::Vector{Float32}, q::Vector{Float32}) = gwstd(x, q)
wstd(x::Vector{Float64}, q::Vector{Float64}) = gwstd(x, q)
wstd(x::Vector{Number}, q::Vector{Number}) = gwstd(x, q)


"""
	mean_std(x, xmin, xmax)
	Returns mean and std for a vector x in the interval between xmin and xmax
"""
function gmean_std(x, xmin, xmax)
    xx = in_range(x, xmin, xmax)
    xm = mean(xx)
    xs = std(xx)
    return xm, xs
end
mean_std(x::Vector{Float32}, xmin::Float32, xmax::Float32) = gmean_std(x, xmin, xmax)
mean_std(x::Vector{Float64}, xmin::Float64, xmax::Float64) = gmean_std(x, xmin, xmax)
mean_std(x::Vector{Int64}, xmin::Int64, xmax::Int64) = gmean_std(x, xmin, xmax)
