using DataFrames

"""
    setunits(df::DataFrame)

Return a dataframe for lor computation where quantities are defined with their units
"""
function setunits(df::DataFrame)
    DataFrame(
        t1 = df.t1 * ns,      # true
        t2 = df.t2 * ns,
        tr1 = df.tr1 * ns,    # reco (smeared take minimum)
        tr2 = df.tr2 * ns,
        ta1 = df.ta1 * ns,    # average (over 5 mimimum in t)
        ta2 = df.ta2 * ns,
        x1 = df.x1 * mm,      # reco position
        x2 = df.x2 * mm,
        y1 = df.y1 * mm,
        y2 = df.y2 * mm,
        z1 = df.z1 * mm,
        z2 = df.z2 * mm,
        xt1 = df.xt1 * mm,      # true position
        xt2 = df.xt2 * mm,
        yt1 = df.yt1 * mm,
        yt2 = df.yt2 * mm,
        zt1 = df.zt1 * mm,
        zt2 = df.zt2 * mm,
        xb1 = df.xb1 * mm,      # position silicon with time stamp x
        xb2 = df.xb2 * mm,
        yb1 = df.yb1 * mm,
        yb2 = df.yb2 * mm,
        zb1 = df.zb1 * mm,
        zb2 = df.zb2 * mm,
    )
end

"""
    deltatime(df::DataFrame, t::String="first")

Return t2 - t1 where t1 and t2 are:
true time stamp (MC), if t : nominal
minimum time stamp of all SiPMs if t : minimum
average time stamp of 5 minimum SiPMs if t: average
"""
function deltatime(df::DataFrame, t::String="first")
    if t     == "minimum"
        return  df.tr1 - df.tr2
    elseif t == "average"
        return df.ta1 - df.ta2
    else
        return df.t1 - df.t2
    end
end

"""
    cdoi(df::DataFrame, position::String="true", nlxe::Number=1.6)

Return the time correction associated to DOI (depth of interaction)
position:
  reco: if using reconstructed values
  true: if using true values

"""
function cdoi(df::DataFrame, position::String="true", nlxe::Number=1.6)
    clxe = SpeedOfLightInVacuum/nlxe

    if position == "reco"
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
    return (dxrb1 - dxrb2)/clxe
end

"""
    dftolor(df::DataFrame, t::String="first", position::String="true", nlxe::Number=1.6)

Take dataframe df and return a vector of MlemLor.
Since MlemLor will be written to file, remove units (e.g, use implicit units, in this
case mm) and transform to Float32
"""
function dftolor(df::DataFrame, t::String="first", position::String="true", nlxe::Number=1.6)

    function tof32(l)
        Float32.(l/mm)
    end

    dt12  = deltatime(df,t)
    dtdoi = cdoi(df,position)
    # compute dx from time and speed of light, ensure that the result is in mm
    dx    = uconvert.(mm, (dt12 - dtdoi) * SpeedOfLightInVacuum)

    if position == "true"
        x1, x2, y1, y2, z1, z2 = df.xt1, df.xt2, df.yt1, df.yt2, df.zt1, df.zt2
    else
        x1, x2, y1, y2, z1, z2 = df.x1,  df.x2,  df.y1,  df.y2,  df.z1,  df.z2
    end

    MlemLor.(tof32(-dx),tof32(x1),tof32(y1),tof32(z1), tof32(x2), tof32(y2), tof32(z2))
end

"""
    ct12(df::DataFrame, t::String="nominal")

Compute raw time difference in ps form analysis dataframe info
- nominal: Use nominal (unsmeared) time recorded by SiPMs.
- minimum: Use the fastest, smeared time, recorded by SiPMs.
- average: Use the average of the fastest time by 5 SiPMs.
"""
function ct12(df::DataFrame, t::String="nominal")
    if t == "minimum"
        return 1000.0*(df.tr1 - df.tr2)
    elseif t == "average"
        return 1000.0*(df.ta1 - df.ta2)
    else
        return 1000.0 * (df.t1 - df.t2)
    end
end

"""
    ctsr(df::DataFrame, reco="nominal")

Return the time difference (in ps) associated to the position of the source (xs).
Notice that this correction is valid to compute the nominal CRT of the system (e.g,
in NEMA3 where sources are at different locations) but should not be used when
estimating DT for a realistic system (e.g, a phantom)

- oreco: Optimized reconstruction. If true, use best (optimized) reconstruction.

"""
function ctsr(df::DataFrame, oreco::Bool=true)
    if oreco == false
        tsr1 = [JPetalo.dxyz([df.xr1[i], df.yr1[i], df.zr1[i]],
                             [df.xs[i], df.ys[i], df.zs[i]]) for i in 1:nrow(df)]/cc
        tsr2 = [JPetalo.dxyz([df.xr2[i], df.yr2[i], df.zr2[i]],
                             [df.xs[i], df.ys[i], df.zs[i]]) for i in 1:nrow(df)]/cc
    else
        tsr1 = [JPetalo.dxyz([df.x1[i], df.y1[i], df.z1[i]],
                             [df.xs[i], df.ys[i], df.zs[i]]) for i in 1:nrow(df)]/cc
        tsr2 = [JPetalo.dxyz([df.x2[i], df.y2[i], df.z2[i]],
                             [df.xs[i], df.ys[i], df.zs[i]]) for i in 1:nrow(df)]/cc
    end
    return 1000.0*(tsr1 - tsr2)
end


"""
    crb(df::DataFrame, nominal::Bool=true, nlxe::Number=1.6)

Return the correction of time of flight associated to the distance travelled by optical
photons in Lxe, from impact point to sensors.

- oreco: Optimized reconstruction. If true, use best (optimized) reconstruction.

"""
function crb(df::DataFrame, oreco::Bool=true, nlxe::Number=1.6)
    clxe = cc/nlxe
    if oreco == false
        trb1 = [JPetalo.dxyz([df.xr1[i], df.yr1[i], df.zr1[i]],
                             [df.xb1[i], df.yb1[i], df.zb1[i]]) for i in 1:nrow(df)]/clxe;
        trb2 = [JPetalo.dxyz([df.xr2[i], df.yr2[i], df.zr2[i]],
                             [df.xb2[i], df.yb2[i], df.zb2[i]]) for i in 1:nrow(df)]/clxe;
    else
        trb1 = [JPetalo.dxyz([df.x1[i], df.y1[i], df.z1[i]],
                             [df.xb1[i], df.yb1[i], df.zb1[i]]) for i in 1:nrow(df)]/clxe;
        trb2 = [JPetalo.dxyz([df.x2[i], df.y2[i], df.z2[i]],
                             [df.xb2[i], df.yb2[i], df.zb2[i]]) for i in 1:nrow(df)]/clxe;
     end
    return 1000.0*(trb1 - trb2)
end

"""
    crt(df::DataFrame, t="nominal", oreco=true, tsr=false)

Return CRT
Value of parameters:

t
- nominal: Use nominal (unsmeared) time recorded by SiPMs.
- minimum: Use the fastest, smeared time, recorded by SiPMs.
- average: Use the average of the fastest time by 5 SiPMs.

- oreco: Optimized reconstruction. If true, use best (optimized) reconstruction.

-tsr : If true, correct for source position.
"""
function crt(df::DataFrame, t="nominal", oreco=true, tsr=false)
    t12 = ct12(df, t)

    tsr12 = ctsr(df, oreco)
    trb12 = crb(df, oreco)

    if tsr
        dt = t12 - tsr12 - trb12
    else
        dt = t12  - trb12
    end

    return t12, dt
end
