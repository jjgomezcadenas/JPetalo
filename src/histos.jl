using Plots
using Statistics
using StatsBase
# histograms

bigp=1.0e+10
bign=-1.0e+10
"""
    digitize(x, bins)

Return the index of the element in the array
"""
function digitize(x::Vector{T}, bins::LinRange{T}) where T
    return searchsortedlast.(Ref(bins), x)
end
#digitize(x::Vector{Float64}, bins::LinRange{Float64}) = gdigitize(x, bins)
#digitize(x::Vector{Float32}, bins::LinRange{Float32}) = gdigitize(x, bins)
#digitize(x::Vector{Real}, bins::LinRange{Number}) = gdigitize(x, bins)


"""
    hist1d(x::Vector{T}, nbins::Integer, xmin::T=bign, xmax::T=bigp)
    hist1d(x::Vector{T}, xs::String, nbins::Integer,
                    xmin::T=bign, xmax::T=bigp; datap = true, fwl=true)
    hist1d(h::Histogram, xs::String; datap = true, markersize=3, fwl=false)

return a 1d histogram and its corresponding graphics (plots)
"""
function hist1d(x::Vector{T}, nbins::Integer, xmin::T=bign, xmax::T=bigp, norm=false) where T
    xx = in_range(x, xmin, xmax)
    dx = (xmax - xmin) / nbins
    bins =[xmin + i * dx for i in 0:nbins]
    h = fit(Histogram, xx, bins)
    if norm
        h = StatsBase.normalize(h, mode=:density)
    end
    return h
end

function hist1d(x::Vector{T}, xs::String, nbins::Integer,
                xmin::T=bign, xmax::T=bigp;
                norm=false, datap = true, markersize=3, fwl=false) where T

    return hist1d(hist1d(x, nbins, xmin, xmax, norm), xs,
                         datap=datap, markersize=markersize, fwl=fwl)
end

function hist1d(h::Histogram, xs::String; datap=true, markersize=3, fwl=false)

    if datap
        yg = h.weights * 1.0
        xg = centers(h)
        p = scatter(xg,yg, yerr = sqrt.(yg), markersize=markersize, legend=false)
        if fwl
            p = plot!(p, xg,yg, yerr = sqrt.(yg), linewidth=1, legend=false)
        end
    else
        p = plot(h, xlabel=xl, yl="frequency")
    end
    xlabel!(xs)
    ylabel!("frequency")

    return h, p
end

function hist1d(h1::Histogram, h2::Histogram, xs::String; markersize=2, norm=false)

    if norm
        h1 = StatsBase.normalize(h1, mode=:density)
        h2 = StatsBase.normalize(h2, mode=:density)
    end

    yg1 = h1.weights * 1.0
    xg1 = centers(h1)
    yg2 = h2.weights * 1.0
    xg2 = centers(h2)

    p1 = scatter(xg1,yg1, yerr = sqrt.(yg1), markersize=markersize)
    p  = scatter!(p1,xg2,yg2, yerr = sqrt.(yg2), markersize=markersize)
    xlabel!(xs)
    ylabel!("frequency")

    return p
end


# function ghist1d(x, nbins, xmin=bign, xmax=bigp)
#     xx = in_range(x, xmin, xmax)
#     h = fit(Histogram, xx, nbins=nbins)
#     return h
# end
# function g2hist1d(x, xl, nbins, xmin=bign, xmax=bigp)
#     h = ghist1d(x, nbins, xmin, xmax)
#     p = plot(h, xlabel=xl, yl="frequency")
#     return h, p
# end
# hist1d(x::Vector{Float64}, nbins::Integer,
#        xmin::Float64=bign, xmax::Float64=bigp) = ghist1d(x, nbins, xmin, xmax)
# hist1d(x::Vector{Float64}, xl::String, nbins::Integer,
#       xmin::Float64=bign, xmax::Float64=bigp) = g2hist1d(x, xl, nbins, xmin, xmax)
#
# hist1d(x::Vector{Float32}, nbins::Integer,
#        xmin::Float32=Float32(bign), xmax::Float32=Float32(bigp)) = ghist1d(x, nbins, xmin, xmax)
# hist1d(x::Vector{Float32}, xl::String, nbins::Integer,
#       xmin::Float32=Float32(bign), xmax::Float32=Float32(bigp)) = g2hist1d(x, xl, nbins, xmin, xmax)
#
# hist1d(x::Vector{Float32}, nbins::Integer,
#        xmin::Float64=bign, xmax::Float64=bigp) = ghist1d(x, nbins, xmin, xmax)
# hist1d(x::Vector{Float32}, xl::String, nbins::Integer,
#       xmin::Float64=bign, xmax::Float64=bigp) = g2hist1d(x, xl, nbins, xmin, xmax)
#
# hist1d(x::Vector{Int64}, nbins::Int64,
#        xmin::Int64=-1000000000, xmax::Int64=1000000000) = ghist1d(x, nbins, xmin, xmax)
# hist1d(x::Vector{Int64}, xl::String, nbins::Int64,
#       xmin::Int64=-1000000000, xmax::Int64=1000000000) = g2hist1d(x, xl, nbins, xmin, xmax)


"""
    hist2d(x,y, nbins, xl, yl)

return a 2d histogram and its corresponding graphics (plots)
"""
function hist2d(x::Vector{T},y::Vector{T}, nbins::Integer,
                xl::String, yl::String,
                xmin::T=bign, xmax::T=bigp,ymin::T=bign, ymax::T=bigp) where T
    function xy(i)
        return diff(h.edges[i])/2 .+ h.edges[i][1:end-1]
    end

    df = DataFrame(x=y,y=x)
    df1 = select_by_column_value_interval(df, "y", xmin, xmax)
    df2 = select_by_column_value_interval(df1, "x", ymin, ymax)
    data = (df2.x, df2.y)
    h = fit(Histogram, data, nbins=nbins)
    ye = xy(1)
    xe = xy(2)
    hm = heatmap(xe, ye, h.weights)
    xlabel!(xl)
    ylabel!(yl)
    return h,hm
end

# function ghist2d(x,y, nbins, xl, yl,xmin=-1e+9, xmax=1e+9,ymin=-1e+9, ymax=1e+9)
#     function xy(i)
#         return diff(h.edges[i])/2 .+ h.edges[i][1:end-1]
#     end
#     df = DataFrame(x=y,y=x)
#     df1 = JPetalo.select_by_column_value_interval(df, "y", xmin, xmax)
#     df2 = JPetalo.select_by_column_value_interval(df1, "x", ymin, ymax)
#     data = (df2.x, df2.y)
#     h = fit(Histogram, data, nbins=nbins)
#     ye = xy(1)
#     xe = xy(2)
#     hm = heatmap(xe, ye, h.weights)
#     xlabel!(xl)
#     ylabel!(yl)
#     return xe,ye,h,hm
# end
#
# hist2d(x::Vector{Float64}, y::Vector{Float64}, nbins::Integer,
#        xl::String, yl::String,
#        xmin::Float64=-1e+9, xmax::Float64=1e+9,
#        ymin::Float64=-1e+9, ymax::Float64=1e+9) = ghist2d(x,y, nbins, xl, yl, xmin, xmax, ymin, ymax)
# hist2d(x::Vector{Float32}, y::Vector{Float32}, nbins::Integer,
#        xl::String, yl::String,
#        xmin::Float32=-1e+9, xmax::Float32=1e+9,
#        ymin::Float32=-1e+9, ymax::Float32=1e+9) = ghist2d(x,y, nbins, xl, yl, xmin, xmax, ymin, ymax)
# hist2d(x::Vector{Int64}, y::Vector{Int64}, nbins::Integer,
#        xl::String, yl::String,
#        xmin::Int64=-1000000000, xmax::Int64=1000000000,
#        ymin::Int64=-1000000000, ymax::Int64=-1000000000) = ghist2d(x,y, nbins, xl, yl, xmin, xmax, ymin, ymax)


"""
    p1df(x, y, nbins)

return a profile DataFrame. This is a DF in which the variable y is
histogrammed as a function of the average of variable x in each bin
"""
function gp1df(x, y, nbins)
    df = DataFrame(x =  x, y = y)                      # create the DF
    bins = LinRange(minimum(x), maximum(x), nbins)     # bins in x
    df[!, "bin"] = digitize(x, bins)
    bin_centers =[0.5 * (bins[i] + bins[i+1]) for i in 1:length(bins) -1]
    bin_width = [bins[i+1] - bins[i] for i in 1:length(bins) -1]
    # mean and std of y in binned x
    ymean = combine(groupby(df, :bin), :y => mean)
    ystd  = combine(groupby(df, :bin), :y => std)
    ndf   = DataFrame(y_mean = ymean[1:end-1, "y_mean"],y_std = ystd[1:end-1, "y_std"])
    ndf[!, "x_mean"] = bin_centers
    ndf[!, "x_std"] = bin_width ./2.0
    return ndf
end

p1df(x::Vector{Float64}, y::Vector{Float64},
     nbins::Integer)  = gp1df(x, y, nbins)
p1df(x::Vector{Float32}, y::Vector{Float32},
     nbins::Integer)  = gp1df(x, y, nbins)
p1df(x::Vector{Number}, y::Vector{Number},
     nbins::Integer)  = gp1df(x, y, nbins)

"""
    centers(h::Histogram)

centers of the histogram
"""
 function centers(h::Histogram)
     edges = collect(h.edges[1])
     return [0.5 *(edges[i] + edges[i+1]) for i in 1:length(edges)-1]
 end
