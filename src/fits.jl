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

function fit_pol2(x::Vector{T},y::Vector{T}) where T
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

function gausg(μ::T, σ::T, C::T) where T
	function gausx(x)
		return C * pdf(Normal(μ, σ,), x)
	end
	return gausx
end

function gausg2(μ1::T, σ1::T, C1::T, μ2::T, σ2::T, C2::T) where T
	function gausx(x)
		return C1 * pdf(Normal(μ1, σ1,), x) + C2 * pdf(Normal(μ2, σ2,), x)
	end
	return gausx
end

@. gauss1fm(x, p) = p[1]* pdf(Normal(p[2], p[3]), x)
@. gauss1(x, p) = p[1]* pdf(Normal(p[2], p[3]), x)
@. gauss2(x, p) = p[1]* pdf(Normal(p[2], p[3]), x) +  p[4]* pdf(Normal(p[5], p[6]), x)

function gaussfm(mu::T) where T
	function gauss(x::Vector{T}, p::Vector{T}) where T
		return @. p[1]* pdf(Normal(mu, p[2]), x)
	end
	return gauss
end

function gauss2fm(mu::T) where T
	function gauss2(x::Vector{T}, p::Vector{T}) where T
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

function cfit(ffit::Function, x::Vector{Float64}, y::Vector{Float64},
	                         p0::Vector{Float64}, lb::Vector{Float64}, ub::Vector{Float64})
	fq = curve_fit(ffit, x, y, p0, lower=lb, upper=ub)
    cfq = coef(fq)
    @debug "coef(fq)" cfq
    return cfq
end

function fit_gauss(h::Histogram)
	c = centers(h)
	w = h.weights * 1.0
	@debug "histo"  w c
	mu, sigma = mean_and_std(c, Weights(h.weights); corrected = false)
	@debug "mu, std" mu, sigma

	# fit parameters lb, ub, po are lower, upper bounds and pars
    lb = [0., mu - 10.0*sigma, sigma/10.0]
    ub = [100.0*sum(w), mu + 10.0*sigma, 10.0*sigma]
    p0_bounds = [sum(w), mu, sigma]
	CC, μ, σ  = cfit(gauss1, c, w, p0_bounds, lb, ub)
	gx = gausg(μ, σ, CC)
	return FGauss([μ], [σ], [CC], h, c, gx.(c), [gx])

end


fit_gauss(x::Vector{Float64}, xmin::Float64, xmax::Float64;
	      bins::Integer=50, norm=false) =fit_gauss(hist1d(x, bins, xmin, xmax, norm))


"""
	fit_gauss_fm(y::Vector{Float64}, xmin::Float64, xmax::Float64, bins=50, fm=0.0)

Fit a gaussian with a fixed mean fm
"""
function fit_gauss_fm(y::Vector{Float64}, xmin::Float64, xmax::Float64;
	                  bins=50, norm=false, fm=0.0)

	# fit the unbinned distribution
	x  =  in_range(y, xmin, xmax)
	σ = std(x)
	@debug "gfit_gauss_fm: σ = $σ"

	# bin distribution
    h = hist1d(x, bins, xmin, xmax, norm)
    c = centers(h)
    w = h.weights *1.0
    @debug "histo w and c"  w c

	# fit parameters lb, ub, po are lower, upper bounds and pars
    lb = [0.0, σ/10.0]
    ub = [1.0e+6*sum(w),  σ * 10.0]
    p0_bounds = [0.0, σ]

	g1 = gaussfm(fm)
	CC, sigma  = cfit(g1, c, w, p0_bounds, lb, ub)

	@debug "CC,  sigma"  CC  sigma
	gx = gausg(fm, sigma, CC)

	return FGauss([fm], [sigma], [CC], h, c, gx.(c), [gx])
end


"""
	fitg1(x, xs, xmin, xmax, xgmin, xgmax; bins=100)

returns the fit and the plot
"""
function fitg1(x, xs, bins, xmin, xmax;
	           xgmin, xgmax, fbins=100, norm=true, fm=0.0)

	h, p = hist1d(x, xs, bins, xmin, xmax, norm=norm)
    fg = fit_gauss_fm(x, xgmin, xgmax, bins=fbins, norm=norm, fm=fm)
	gx = fg.g[1]
	X = centers(h)
    p = plot!(p, X, gx.(X), lw=2, legend=false)
	xlabel!(xs)
	return fg, p
end

"""
	gfit_gauss2_cmean(y, xmin, xmax, bins, sigmas, cs, cmean=0.0)

Fit a double gaussian (with sigmas -->[sigma1, sigma2] cs -->[c1, c2] )
and a cmean to data.
"""
function gfit_gauss2_cmean(y::Vector{Float64}, xmin::Float64, xmax::Float64,
	                       bins::Integer, sigmas::Vector{Float64}, cs::Vector{Float64},
						   norm=false, cmean=0.0)

    x =  in_range(y, xmin, xmax)
    h = hist1d(x, bins, xmin, xmax, norm)
    c = centers(h)
    w = h.weights
    @debug "histo centers and weights in full region"  w c

    g2 = gauss2fm(cmean)
    # fit parameters lb, ub, po are lower, upper bounds and pars

    lb = [cs[1]/100.0, sigmas[1]/5.0, cs[2]/100.0, sigmas[2]/5.0]
    ub = [cs[1]*100.0, sigmas[1]*5.0, cs[2]*100.0, sigmas[2]*5.0]
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
    #return (sigma1 = sigma1, sigma2 = sigma2, C1 = C1, C2=C2,
    #        h = h, xg = c, yg = gx.(c), gx = gx, gx1=gx1, gx2=gx2)

	return FGauss([cmean, cmean], [sigma1, sigma2], [C1, C2],
			       h, c, gx.(c), [gx, gx1, gx2])
    end

"""
	fit_2gauss_cmean(data, gp, g1p, g2p, cm)

Fit two gaussian with common mean
"""
function fit_2gauss_cmean(data, gp, g1p, g2p, cm, norm=false)
    gf1 = fit_gauss_fm(data, g1p.xmin,g1p.xmax, bins=g1p.nbin, norm=norm, fm=cm)
    @debug gf1
    gf2 = fit_gauss_fm(data, g2p.xmin,g2p.xmax,bins=g2p.nbin, norm= norm, fm=cm)
    @debug gf2
    gf = gfit_gauss2_cmean(data, gp.xmin,gp.xmax,gp.nbin,
	                       [gf1.std[1], gf2.std[1]], [gf1.C[1], gf2.C[1]], norm)
    @debug gf
    return gf
end


# fit_gauss2(x::Vector{Float64},
#            xmin::Vector{Float64},
# 		   xmax::Vector{Float64},
# 		   bins::Vector{Int64}) = gfit_gauss2(x, xmin, xmax, bins)


"""
    fitg2(x, xs, xmin, xmax, xg1min, xg1max, xg2min, xg2max, xgmin, xgmax; bins=100)
Fits 2 gaussians with common mean (0 by default) to vector x.
"""
function fitg2(x, xs, bins, xmin, xmax;
	           xg1min, xg1max, xg2min, xg2max, xgmin, xgmax, cm=0.0,
         	   g1bins=100, g2bins=100, gbins=100, norm=true)

	hp, p = hist1d(x, xs, bins, xmin, xmax, norm=norm)

    g1p = (xmin = xg1min, xmax = xg1max, nbin=g1bins)
    g2p = (xmin= xg2min, xmax=  xg2max, nbin=g2bins)
    gp  = (xmin= xgmin, xmax=  xgmax, nbin=gbins)

    fg  = fit_2gauss_cmean(x, gp, g1p, g2p, cm, norm)
	gx  = fg.g[1]
	gx1 = fg.g[2]
	gx2 = fg.g[3]

    p = plot!(p,fg.X, fg.Y, lw=2)
    p = plot!(p,fg.X, gx1.(fg.X), lw=1)
    p = plot!(p,fg.X, gx2.(fg.X), lw=1)
    #xlabel!(xs)
    return fg, p
end


"""
    fit_profile(df1, c1, c2, tx1, ty1, fit="pol1", bins=25)
Create and fit a profile with pol1 or poli2 functions.
Return fit parameters, fit function and plot
"""
function fit_profile(df1, c1, c2, tx1, ty1, fit="pol1", bins=25)
    function gf1(ct, fit)
        function f1(z)
            return ct[1] + ct[2] * z
        end

        function f2(z)
            return ct[1] + ct[2] * z + ct[3] * z^2
        end

        if fit == "pol1"
            return f1
        end
        return f2
    end

    pdf1 = p1df(df1[!,c1],df1[!,c2], bins)

    if fit == "pol1"
        lt1, pt1, ct1 = lfit(pdf1)
    else
        ct1 = fit_pol2(pdf1.x_mean, pdf1.y_mean)
    end

    ff1 = gf1(ct1, fit)

    p1 = plot(pdf1.x_mean,pdf1.y_mean, yerror=pdf1.y_std, shape = :circle, color = :black, legend=false)
    p1 = plot!(p1, pdf1.x_mean, ff1.(pdf1.x_mean))
    xlabel!(tx1)
    ylabel!(ty1)

   return ct1, ff1, p1
end
