using DataFrames
using Glob
using CSV
using DrWatson
#using LsqFit
#using GLM

"""
    readdf(dir)

Reads all csv files found in dir and returns a single DataFrame.
"""
function readdf(dir)
    drx = datadir(dir)
    files = glob("*.csv",drx)
    dfs =[DataFrame(CSV.File(file)) for file in files]
    evtdf=vcat(dfs...)
end

"""
    selectinterval(df, column1, column2, xmin, xmax)

Select interval over 2 hemisphere columns (e.g, q1, q2) for df
"""
function selectinterval(df, column1, column2, xmin, xmax)
    e1 = JPetalo.select_by_column_value_interval(df, column1, xmin,xmax)
    JPetalo.select_by_column_value_interval(e1, column2, xmin,xmax)
end


"""
    plotreso(dfx, r1x, tx1, ty1, xmin, xmax, bins=150)
Plots the resolution of estimator r1x wrt true radius r1

"""
function plotreso(dfx, r1x, tx1, ty1, xmin, xmax, bins=150)
    h1,p1 = hist2d(dfx.r1,r1x, bins, tx1, ty1)
    h2,p2       = hist1d(dfx.r1 - r1x, "DX", bins, xmin, xmax)
    plot(p1, p2,  layout= (1, 2), legend=false, fmt = :png,  size = (1000, 400),
        left_margin=5Plots.mm, right_margin=1Plots.mm, bottom_margin=5Plots.mm)
end


"""
    ploth2d(df1, c1, c2, tx1, ty1, xmin1, xmax1, xmin2, xmax2, bins=150)

Return a 2d histogram of variables of columns c1 and c2 of df1.
"""
function ploth2d(df1, c1, c2, tx1, ty1, xmin1, xmax1, xmin2, xmax2, bins=150)
    h1,p1 = hist2d(df1[!,c1],df1[!,c2], bins, tx1, ty1, xmin1, xmax1, xmin2, xmax2)
    plot(p1, layout= (1, 1), legend=false, fmt = :png,
        left_margin=5Plots.mm, right_margin=1Plots.mm, bottom_margin=5Plots.mm)
end


"""
    profile(df1, c1, c2, tx1, ty1, bins=25)
Return a profile histogram plot
"""
function profile(df1, c1, c2, tx1, ty1, bins=25)
    pdf1 = p1df(df1[!,c1],df1[!,c2], bins)
    p1 = plot(pdf1.x_mean,pdf1.y_mean, yerror=pdf1.y_std, shape = :circle, color = :black, legend=false)
    xlabel!(tx1)
    ylabel!(ty1)
    plot(p1, layout= (1, 1), legend=false, fmt = :png,
        left_margin=5Plots.mm, right_margin=1Plots.mm, bottom_margin=5Plots.mm)
end
