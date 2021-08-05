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


function writemdf(dir, file, df)
    path = string(dir,"/", file)
    CSV.write(path, df)
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
function plotreso(r1t, r1x, tx1, ty1, xs, xmin, xmax, nbins=150)
    h1,p1 = hist2d(r1x, r1t, nbins, tx1, ty1)
    h2,p2 = hist1d(r1t - r1x, xs, nbins, xmin, xmax)
    plot(p1, p2,  layout= (1, 2), legend=false, fmt = :png,  size = (1000, 400),
        left_margin=5Plots.mm, right_margin=1Plots.mm, bottom_margin=5Plots.mm)
end

function saveplot(p, filename)
    fout = string(datadir("plots"),"/", filename)
    png(p, fout)
end


function plot_and_save(p1,p2, tit="", filename="")
    if tit != ""
        p = plot(p1, p2, layout= (1, 2), title=tit, legend=false, fmt = :png, size = (1000, 400),
            left_margin=5Plots.mm, right_margin=1Plots.mm, bottom_margin=5Plots.mm)
    else
        p = plot(p1, p2, layout= (1, 2), legend=false, fmt = :png, size = (1000, 400),
            left_margin=5Plots.mm, right_margin=1Plots.mm, bottom_margin=5Plots.mm)
    end
    if filename != ""
        fout = string(datadir("plots"),"/", filename)
        png(fout)
    end
    return p
end

function q1vsq2(df; tit="", filename="")
    h1,p1 = JPetalo.hist2d(df.q1, df.q2, 150, "q1 (pes)","q2 (pes)", 100., 3500.,100., 3500.)
    h2,p2 = JPetalo.hist1d(df.q1, "q1 (pes)", 100, 100.0, 3500.);
    plot(p1, p2, layout= (1, 2), legend=false, fmt = :png, size = (1000, 400),
        left_margin=5Plots.mm, right_margin=1Plots.mm, bottom_margin=5Plots.mm)

    plot_and_save(p1,p2, tit, filename)
end


function r1q1(df, qmin=1500.0, qmax=3000.0; tit="", filename="")
    h1,p1 = JPetalo.hist2d(df.q1, df.r1, 150, "q1 (pes) ","r1 (mm)",qmin, qmax, 350., 450.)
    h2,p2 = JPetalo.hist1d(df.q1, "q1", 100, qmin, qmax,)
    plot(p1, p2, layout= (1, 2), legend=false, fmt = :png, size = (1000, 400),
        left_margin=5Plots.mm, right_margin=1Plots.mm, bottom_margin=5Plots.mm)
    plot_and_save(p1,p2, tit, filename)
end


function zstd(df, zmin=0.0, zmax=40.0; tit="", filename="")
    h1,p1 = JPetalo.hist2d(df.zstd1, df.r1, 150, "σz (mm) ","1 (mm)",zmin, zmax,350., 450.)
    h2,p2 = JPetalo.hist1d(df.zstd1, "σz (mm)", 100, zmin, zmax)
    plot(p1, p2, layout= (1, 2), legend=false, fmt = :png, size = (1000, 400),
        left_margin=5Plots.mm, right_margin=1Plots.mm, bottom_margin=5Plots.mm)
    plot_and_save(p1,p2, tit, filename)
end


function phistd(df, phimin=0.0, phimax=0.1; tit="", filename="")
    h1,p1 = JPetalo.hist2d(df.phistd1, df.r1, 150, "σϕ (mm) ","r1 (mm)",phimin, phimax,350., 450.)
    h2,p2 = JPetalo.hist1d(df.phistd1, "σϕ (mm)", 100, phimin, phimax)
    plot(p1, p2, layout= (1, 2), legend=false, fmt = :png, size = (1000, 400),
        left_margin=5Plots.mm, right_margin=1Plots.mm, bottom_margin=5Plots.mm)

    plot_and_save(p1,p2, tit, filename)
end


function nfit_profile(df,cx,cy,tx,ty,pol)
    cq, fq, pq = JPetalo.fit_profile(df, cx,cy,tx,ty,pol)
    println("for ",cx,"=f(",cy,"): fit parameters =",cq )
    return cq, fq, pq
end


function nplot_profile(p; filename="")
    p= plot(p,  layout= (1, 1), legend=false, fmt = :png,
        left_margin=5Plots.mm, right_margin=1Plots.mm, bottom_margin=5Plots.mm)
    if filename != ""
        fout = string(datadir("plots"),"/", filename)
        png(p,fout)
    end
    p
end
