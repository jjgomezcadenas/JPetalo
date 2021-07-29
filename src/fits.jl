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

@. gauss1fm(x, p) = p[1]* pdf(Normal(p[2], p[3]), x)
@. gauss1(x, p) = p[1]* pdf(Normal(p[2], p[3]), x)
@. gauss2(x, p) = p[1]* pdf(Normal(p[2], p[3]), x) +  p[4]* pdf(Normal(p[5], p[6]), x)

function gaussfm(mu)
	function gauss(x,p)
		return @. p[1]* pdf(Normal(mu, p[2]), x)
	end
	return gauss
end

function gauss2fm(mu)
	function gauss2(x,p)
		return @. p[1]* pdf(Normal(mu, p[2]), x) + p[3]* pdf(Normal(mu, p[4]), x)
	end
	return gauss2
end

function gauss1fm(mu)
    function gauss1(x,p)
        return @. p[1]* pdf(Normal(mu, p[2]), x)
    end
    return gauss1
end

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

function gfit_gauss_fm(y, xmin, xmax, bins=25, fm=0.0)

	# fit the unbinned distribution
	x  =  in_range(y, xmin, xmax)
	σ = std(x)
	@debug "gfit_gauss_fm: σ = $σ"

	# bin distribution
    h = hist1d(x, bins, xmin, xmax)
    c = centers(h)
    w = h.weights
    @debug "histo w and c"  w c

	# fit parameters lb, ub, po are lower, upper bounds and pars
    lb = [sum(w)/50.0, σ/10.0]
    ub = [sum(w),  σ * 10.0]
    p0_bounds = [sum(w)/10.0, σ]

	g1 = gaussfm(fm)
	CC, sigma  = cfit(g1, c, w, p0_bounds, lb, ub)

	@debug "CC,  sigma"  CC  sigma
	gx = gausg(fm, sigma, CC)

	return FGauss([fm], [sigma], [CC], h, c, gx.(c), [gx])
end


function gfit_gauss(y, xmin, xmax, bins=25, ff=gauss1)

	# fit the unbinned distribution
	x  =  in_range(y, xmin, xmax)
	ft =  fit(Normal,x)
	@debug "initial mean and std" ft.μ ft.σ

	# bin distribution
    h = hist1d(x, bins, xmin, xmax)
    c = centers(h)
    w = h.weights
    @debug "histo w and c"  w c

	# fit parameters lb, ub, po are lower, upper bounds and pars
    lb = [sum(w)/50.0, ft.μ - 5*ft.σ, ft.σ/10.0]
    ub = [sum(w), ft.μ + 5*ft.σ, 10*ft.σ]
    p0_bounds = [sum(w)/10.0, ft.μ, ft.σ]
	CC, mu, sigma  = cfit(ff, c, w, p0_bounds, lb, ub)
	gx = gausg(mu, sigma, CC)

	return FGauss([mu], [sigma], [CC], h, c, gx.(c), [gx])
end

"""
	gfit_xgauss(y, xmin, xmax, bins=25, fm=0.0)
	fits a gaussian with a fixed mean
"""
function gfit_xgauss(y, xmin, xmax, bins=25, fm=0.0)

    # fit the unbinned distribution
    x  =  in_range(y, xmin, xmax)
    mu, sigma =  mean_and_std(x)
    @debug "initial std" sigma

    # bin distribution
    h = hist1d(x, bins, xmin, xmax)
    c = centers(h)
    w = h.weights
    @debug "histo w and c"  w c

    # fit parameters lb, ub, po are lower, upper bounds and pars
    lb = [sum(w)/100000.0,  sigma/10.0]
    ub = [sum(w)*100000.0, 10*sigma]
    p0_bounds = [sum(w)/10.0, sigma]

    ff = gauss1fm(fm)
    fq = curve_fit(ff, c, w, p0_bounds, lower=lb, upper=ub)
    cfq = coef(fq)
    @debug "coef(fq)" cfq
    CC, sigma  = coef(fq)

    @debug "coef(fq)" CC  sigma
    gx = gausg(fm, sigma, CC)

    return (mu = fm, sigma = sigma, C = CC,
            h = h, xg = c, yg = gx.(c), gx = gx)
end

"""
	gfit_gauss2_cmean(y, xmin, xmax, bins, sigmas, cs, cmean=0.0)

Fit a double gaussian (with sigmas -->[sigma1, sigma2] cs -->[c1, c2] )
and a cmean to data.
"""
function gfit_gauss2_cmean(y, xmin, xmax, bins, sigmas, cs, cmean=0.0)

    x =  in_range(y, xmin, xmax)
    h = hist1d(x, bins, xmin, xmax)
    c = centers(h)
    w = h.weights
    @debug "histo centers and weights in full region"  w c

    g2 = gauss2fm(cmean)
    # fit parameters lb, ub, po are lower, upper bounds and pars

    lb = [cs[1]/10.0, sigmas[1]/3.0, cs[2]/10.0, sigmas[2]/3.0]
    ub = [cs[1]*10.0, sigmas[1]*3.0, cs[2]*10.0, sigmas[2]*3.0]
    p0_bounds = [cs[1], sigmas[1], cs[2], sigmas[2]]

    @debug "pars" p0 lb ub
    # fit double gaussian
    fq = curve_fit(g2, c, w, p0_bounds, lower=lb, upper=ub)
    C1, sigma1, C2,  sigma2   = coef(fq)
    @debug "C1 sigma1 C2 sigma2" C1 sigma1 C2  sigma2

    #
    gx = gausg2(cmean, sigma1, C1, cmean, sigma2, C2)
    gx1 = gausg(cmean, sigma1, C1)
    gx2 = gausg(cmean, sigma2, C2)
    return (sigma1 = sigma1, sigma2 = sigma2, C1 = C1, C2=C2,
            h = h, xg = c, yg = gx.(c), gx = gx, gx1=gx1, gx2=gx2)
    end

"""
	fit_2gauss_cmean(data, gp, g1p, g2p, cm)

Fit two gaussian with common mean
"""
function fit_2gauss_cmean(data, gp, g1p, g2p, cm)
    gf1 = gfit_xgauss(data, g1p.xmin,g1p.xmax,g1p.nbin, cm)
    @debug gf1
    gf2 = gfit_xgauss(data, g2p.xmin,g2p.xmax,g2p.nbin, cm)
    @debug gf2
    gf = gfit_gauss2_cmean(data, gp.xmin,gp.xmax,gp.nbin, [gf1.sigma, gf2.sigma], [gf1.C, gf2.C])
    @debug gf
    return gf
end

fit_gauss(x::Vector{Float64},
        xmin::Float64, xmax::Float64, bins::Integer=50) = gfit_gauss(x, xmin, xmax, bins)
fit_gauss(x::Vector{Float32},
        xmin::Float32, xmax::Float32, bins::Integer=50) = gfit_gauss(x, xmin, xmax, bins)
fit_gauss(h::Histogram) = hfit_gauss(h::Histogram)
fit_gauss2(x::Vector{Float64},
           xmin::Vector{Float64},
		   xmax::Vector{Float64},
		   bins::Vector{Int64}) = gfit_gauss2(x, xmin, xmax, bins)
