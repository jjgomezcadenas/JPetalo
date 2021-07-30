using DataFrames

@enum Dtsel dtfirst dtminimum dtaverage
@enum Possel postrue posreco

"""
    setunits(df::DataFrame)

Return a dataframe for lor computation where quantities are defined with their units
"""

function setunits(df::DataFrame)
    DataFrame(
        nsipm1 = df.nsipm1,     # number of sipms in cluster
        nsipm2 = df.nsipm2,

        q1 = df.q1,             # q in pes
        q2 = df.q2,

        r1  = df.r1  * mm,      # best radius
        r1q = df.r1q * mm,      # raw radius from q
        r2  = df.r2  * mm,
        r2q = df.r2q * mm,

        t1  = df.t1  * ns,      # time true
        t2  = df.t2  * ns,
        ta1 = df.ta1 * ns,      # average (over 5 mimimum in t)
        ta2 = df.ta2 * ns,
        tr1 = df.tr1 * ns,      # reco (smeared take minimum)
        tr2 = df.tr2 * ns,

        ux = df.ux,             # unit direction vector of gammas
        uy = df.uy,
        uz = df.uz,

        x1  = df.x1  * mm,      # best reco position
        x2  = df.x2  * mm,
        xb1 = df.xb1 * mm,      # position of sipm that gives time stamp
        xb2 = df.xb2 * mm,
        xr1 = df.xr1 * mm,      # reco position
        xr2 = df.xr2 * mm,
        xs  = df.xs  * mm,      # position of source
        xt1 = df.xt1 * mm,      # true position
        xt2 = df.xt2 * mm,

        y1  = df.y1  * mm,
        y2  = df.y2  * mm,
        yb1 = df.yb1 * mm,
        yb2 = df.yb2 * mm,
        yr1 = df.yr1 * mm,
        yr2 = df.yr2 * mm,
        ys  = df.ys  * mm,
        yt1 = df.yt1 * mm,
        yt2 = df.yt2 * mm,

        z1  = df.z1  * mm,
        z2  = df.z2  * mm,
        zb1 = df.zb1 * mm,
        zb2 = df.zb2 * mm,
        zr1 = df.zr1 * mm,
        zr2 = df.zr2 * mm,
        zs  = df.zs  * mm,
        zt1 = df.zt1 * mm,
        zt2 = df.zt2 * mm,

    )
end

"""
    deltatime(df::DataFrame, t::Dtsel=dtfirst)
    @enum Dtsel dtfirst dtminimum dtaverage

Return t2 - t1 where t1 and t2 are:
dtfirst for nominal (true) time
dtminimum for minimum time stamp of all SiPMs
dtaverage for average time stamp of 5 minimum SiPMs
"""
function deltatime(df::DataFrame, t::Dtsel=dtfirst)
    if t     == dtminimum
        return  df.tr2 - df.tr1
    elseif t == dtaverage
        return df.ta2 - df.ta1
    else
        return uconvert.(ps, df.t2 - df.t1)

    end
end

"""
    cdoi(df::DataFrame, position::Possel=postrue, nlxe::Number=1.6)
    @enum Possel postrue posreco

Return the time correction associated to DOI (depth of interaction)
position:
  posreco: if using reconstructed values
  postrue: if using true values

"""
function cdoi(df::DataFrame, position::Possel=postrue, nlxe::Number=1.6)
    clxe = SpeedOfLightInVacuum/nlxe

    if position == posreco
        dxrb1 = [JPetalo.dxyz([df.x1[i], df.y1[i], df.z1[i]],
                              [df.xb1[i], df.yb1[i], df.zb1[i]]) for i in 1:nrow(df)]

        dxrb2 = [JPetalo.dxyz([df.x2[i], df.y2[i], df.z2[i]],
                              [df.xb2[i], df.yb2[i], df.zb2[i]]) for i in 1:nrow(df)]
    else
        dxrb1 = [JPetalo.dxyz([df.xt1[i], df.yt1[i], df.zt1[i]],
                              [df.xb1[i], df.yb1[i], df.zb1[i]]) for i in 1:nrow(df)]

        dxrb2 = [JPetalo.dxyz([df.xt2[i], df.yt2[i], df.zt2[i]],
                              [df.xb2[i], df.yb2[i], df.zb2[i]]) for i in 1:nrow(df)]

    end
    return uconvert.(ps, (dxrb2 - dxrb1)/clxe)
end


function ctsr(df::DataFrame, position::Possel=postrue)
    cc = float(SpeedOfLightInVacuum)
    if position == postrue
        tsr1 = [JPetalo.dxyz([df.xt1[i], df.yt1[i], df.zt1[i]],
                             [df.xs[i], df.ys[i], df.zs[i]]) for i in 1:nrow(df)]/cc
        tsr2 = [JPetalo.dxyz([df.xt2[i], df.yt2[i], df.zt2[i]],
                             [df.xs[i], df.ys[i], df.zs[i]]) for i in 1:nrow(df)]/cc

    else
        tsr1 = [JPetalo.dxyz([df.x1[i], df.y1[i], df.z1[i]],
                             [df.xs[i], df.ys[i], df.zs[i]]) for i in 1:nrow(df)]/cc
        tsr2 = [JPetalo.dxyz([df.x2[i], df.y2[i], df.z2[i]],
                             [df.xs[i], df.ys[i], df.zs[i]]) for i in 1:nrow(df)]/cc
    end
    return uconvert.(ps, tsr2 - tsr1)
end

"""
    crt(dfu, dtsel=JPetalo.dtfirst, posel=JPetalo.postrue)

Return the CRT of the system
"""
function crt(dfu, dtsel=dtfirst, posel=postrue,
            xmin1= -500.0, xmax1= -50.0, nbin1=50,
            xmin2= -50.0,  xmax2=  50.0, nbin2=50,
            xmin= -500.0,  xmax=  500.0, nbin=100)

    dt12 = deltatime(dfu, dtsel)
    t12 = dt12./ps

    dtrb12 = cdoi(dfu, posel)
    trb12 = dtrb12 ./ps

    dtrb12 = cdoi(dfu, posel)
    trb12 = dtrb12 ./ps

    dt = t12 - tsr12 - trb12

    g1p = (xmin= xmin1, xmax= xmax1, nbin=nbin1)
    g2p = (xmin= xmin2, xmax=  xmax2, nbin=nbin2)
    gp  = (xmin= xmin, xmax=  xmax, nbin=nbin)

    fdt = JPetalo.fit_2gauss_cmean(dt, gp, g1p, g2p, 0.0)
    @info "dt: sigma" fdt.sigma1 fdt.sigma2
    return t12, dt, fdt

end

"""
    dftolor(df::DataFrame, t::Dtsel=dtfirst, position::Possel=postrue, nlxe::Number=1.6)

Take dataframe df and return a vector of MlemLor.
Since MlemLor will be written to file, remove units (e.g, use implicit units, in this
case mm) and transform to Float32
"""
function dftolor(df::DataFrame, t::Dtsel=dtfirst, position::Possel=postrue, nlxe::Number=1.6)

    function tof32(l)
        Float32.(l/mm)
    end

    dt12  = deltatime(df,t)
    dtdoi = cdoi(df,position)
    # compute dx from time and speed of light, ensure that the result is in mm
    dx    = uconvert.(mm, (dt12 - dtdoi) * SpeedOfLightInVacuum)

    if position == postrue
        x1, x2, y1, y2, z1, z2 = df.xt1, df.xt2, df.yt1, df.yt2, df.zt1, df.zt2
    else
        x1, x2, y1, y2, z1, z2 = df.x1,  df.x2,  df.y1,  df.y2,  df.z1,  df.z2
    end

    MlemLor.(tof32(dx),tof32(x1),tof32(y1),tof32(z1), tof32(x2), tof32(y2), tof32(z2))
end
