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
    plotreso(r1t, r1x, tx1, ty1, "xs", xmin, xmax, bins=150)
Plots the resolution of estimator r1x wrt true radius r1

"""
function plotreso(r1t, r1x, tx1, ty1, xmin, xmax, bins=150)
    h1,p1 = hist2d(r1x, r1t, bins, tx1, ty1)
    h2,p2 = hist1d(r1t - r1x, xs, bins, xmin, xmax)
    plot(p1, p2,  layout= (1, 2), legend=false, fmt = :png,  size = (1000, 400),
        left_margin=5Plots.mm, right_margin=1Plots.mm, bottom_margin=5Plots.mm)
end
