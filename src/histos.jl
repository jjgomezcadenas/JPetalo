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
function gdigitize(x, bins)
    return searchsortedlast.(Ref(bins), x)
end
digitize(x::Vector{Float64}, bins::LinRange{Float64}) = gdigitize(x, bins)
digitize(x::Vector{Float32}, bins::LinRange{Float32}) = gdigitize(x, bins)
digitize(x::Vector{Real}, bins::LinRange{Number}) = gdigitize(x, bins)


"""
    hist1d(x, nbins, xl)

return a 1d histogram and its corresponding graphics (plots)
"""
function ghist1d(x, nbins, xmin=bign, xmax=bigp)
    xx = in_range(x, xmin, xmax)
    h = fit(Histogram, xx, nbins=nbins)
    return h
end
function g2hist1d(x, xl, nbins, xmin=bign, xmax=bigp)
    h = ghist1d(x, nbins, xmin, xmax)
    p = plot(h, xlabel=xl, yl="frequency")
    return h, p
end

hist1d(x::Vector{Float64}, nbins::Integer,
       xmin::Float64=bign, xmax::Float64=bigp) = ghist1d(x, nbins, xmin, xmax)
hist1d(x::Vector{Float64}, xl::String, nbins::Integer,
      xmin::Float64=bign, xmax::Float64=bigp) = g2hist1d(x, xl, nbins, xmin, xmax)

hist1d(x::Vector{Float32}, nbins::Integer,
       xmin::Float32=Float32(bign), xmax::Float32=Float32(bigp)) = ghist1d(x, nbins, xmin, xmax)
hist1d(x::Vector{Float32}, xl::String, nbins::Integer,
      xmin::Float32=Float32(bign), xmax::Float32=Float32(bigp)) = g2hist1d(x, xl, nbins, xmin, xmax)

hist1d(x::Vector{Float32}, nbins::Integer,
       xmin::Float64=bign, xmax::Float64=bigp) = ghist1d(x, nbins, xmin, xmax)
hist1d(x::Vector{Float32}, xl::String, nbins::Integer,
      xmin::Float64=bign, xmax::Float64=bigp) = g2hist1d(x, xl, nbins, xmin, xmax)

hist1d(x::Vector{Int64}, nbins::Int64,
       xmin::Int64=-1000000000, xmax::Int64=1000000000) = ghist1d(x, nbins, xmin, xmax)
hist1d(x::Vector{Int64}, xl::String, nbins::Int64,
      xmin::Int64=-1000000000, xmax::Int64=1000000000) = g2hist1d(x, xl, nbins, xmin, xmax)


"""
    hist2d(x,y, nbins, xl, yl)

return a 2d histogram and its corresponding graphics (plots)
"""
function ghist2d(x,y, nbins, xl, yl,xmin=-1e+9, xmax=1e+9,ymin=-1e+9, ymax=1e+9)
    function xy(i)
        return diff(h.edges[i])/2 .+ h.edges[i][1:end-1]
    end
    df = DataFrame(x=y,y=x)
    df1 = JPetalo.select_by_column_value_interval(df, "y", xmin, xmax)
    df2 = JPetalo.select_by_column_value_interval(df1, "x", ymin, ymax)
    data = (df2.x, df2.y)
    h = fit(Histogram, data, nbins=nbins)
    ye = xy(1)
    xe = xy(2)
    hm = heatmap(xe, ye, h.weights)
    xlabel!(xl)
    ylabel!(yl)
    return xe,ye,h,hm
end

hist2d(x::Vector{Float64}, y::Vector{Float64}, nbins::Integer,
       xl::String, yl::String,
       xmin::Float64=-1e+9, xmax::Float64=1e+9,
       ymin::Float64=-1e+9, ymax::Float64=1e+9) = ghist2d(x,y, nbins, xl, yl, xmin, xmax, ymin, ymax)
hist2d(x::Vector{Float32}, y::Vector{Float32}, nbins::Integer,
       xl::String, yl::String,
       xmin::Float32=-1e+9, xmax::Float32=1e+9,
       ymin::Float32=-1e+9, ymax::Float32=1e+9) = ghist2d(x,y, nbins, xl, yl, xmin, xmax, ymin, ymax)
hist2d(x::Vector{Int64}, y::Vector{Int64}, nbins::Integer,
       xl::String, yl::String,
       xmin::Int64=-1000000000, xmax::Int64=1000000000,
       ymin::Int64=-1000000000, ymax::Int64=-1000000000) = ghist2d(x,y, nbins, xl, yl, xmin, xmax, ymin, ymax)


"""
    p1df(x, y, nbins)

return a profile DataFrame. This is a DF in which the variable y is
histogrammed as a function of the average of variable x in each bin
"""
# function gp1df(x, y, nbins)
#     df = DataFrame(x =  x, y = y)                      # create the DF
#     bins = LinRange(minimum(x), maximum(x), nbins)     # bins in x
#     # bin x
#     df[!, "bin"] = digitize(x, bins)
#     bin_centers =[0.5 * (bins[i] + bins[i+1]) for i in 1:length(bins) -1]
#     bin_width = [bins[i+1] - bins[i] for i in 1:length(bins) -1]
#     # mean and std of y in binned x
#     ymean = combine(groupby(df, :bin), :y => mean)
#     ystd  = combine(groupby(df, :bin), :y => std)
#     ndf   = innerjoin(ymean, ystd, on = :bin)
#     ndf[!, "x_mean"] = bin_centers
#     ndf[!, "x_std"] = bin_width ./2.0
#     return ndf
# end
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
