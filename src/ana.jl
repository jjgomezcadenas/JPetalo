using DataFrames
using StatsModels
using Clustering
using Statistics

# Selection

"""
    primary_in_lxe(verticesdf::DataFrame)

Select primary photons in LXe.

Row │ event_id track_id parent_id x y z t moved pre_KE post_KE deposited process_id  volume_id
LXe           : volume_id = 0
Primary       : parent_id = 0
"""
function primary_in_lxe(verticesdf::DataFrame)
    vlxe      = select_by_column_value(verticesdf, "volume_id", 0)
    vlxephepr = select_by_column_value(vlxe, "parent_id", 0)
    return vlxephepr
end


"""
	select_truehit(th::TrueHits, index::Integer)

Returns TrueHit corresponding to index from a vector of TrueHits

"""
function select_truehit(th::TrueHits, index::Integer)
	eid = th.event_id[index]
	x   = th.x[index]
	y   = th.y[index]
	z   = th.z[index]
	t   = th.t[index]
	e   = th.e[index]
	return TrueHit(eid,x,y,z,t,e)
end


#Reco
# """
# 	reco_hits(event::Int64, ecut::Float32,
# 			  Qdf::DataFrame, sxyzdf::DataFrame)
#
# Return a RecoHits DataFrame.
# Qdf DF -> pdf.total_charge (event_id, sensor_id, charge)
# sxyzdf ->pdf.sensor_xyz (sensor_id, x, y, z)
# RecoHits DF: (x,y,z,q) for each SiPM in event.
# """

function rhits(event::Integer, ecut::Number, pde::Number,
               Qdf::DataFrame, sxyzdf::DataFrame)

    # select the event
    qdf = select_by_column_value(Qdf, "event_id", event)
    # multiply vector of charges by PDE and add column to the DF
    Q = Float32.(qdf.charge * pde)
    qdf[!,"Q"] = Q
    # Select SiPMs with charge (q x PDE) about ecut
    qdfQ   = qdf[qdf.Q.>ecut,:]
    return sipm_xyzq(qdfQ, sxyzdf)
end


"""
	sipmsel(hdf::DataFrame)

Return two data frames, separating the SiPMs in two groups depending on
the sign of the angle with the SiPM of max charge
"""
function sipmsel(hdf::DataFrame)
	simax = find_xyz_sipm_qmax(hdf)
	npr   = xyz_dot(hdf, simax)
	mask =[n>0 ? true : false for n in npr]
	return hdf[(npr.>0), :], hdf[(npr.<0), :]
end


"""
	ksipmsel(hdf::DataFrame, ka::Vector{Int64})

Return two data frames, separating the SiPMs in two groups depending on
the value of vector ka (1 or 2)
"""
function ksipmsel(hdf::DataFrame, ka::Vector{Int64})
	return hdf[(ka.==2), :], hdf[(ka.==1), :]
end


"""
	baricenter(hdf::DataFrame)
	returns the barycenter of a cluster of hits
"""
function baricenter(hdf::DataFrame)
	function xq(hdf::DataFrame, pos::String)
		return sum(hdf[!,pos] .* hdf.q) / qt
	end
	qt = sum(hdf.q)
	return Hit(xq(hdf, "x"), xq(hdf, "y"), xq(hdf, "z"), qt)
end


# kmeans algorithm
"""
	get_hits_as_matrix(hitdf::DataFrame)

Given the dataframe hitdf (with fields x,y,z), return the
underlying matrix
"""
function get_hits_as_matrix(hitdf::DataFrame)
	f = @formula(0 ~ x + y + z)
	f2 = apply_schema(f, schema(f, hitdf))
	resp, pred = modelcols(f2, hitdf)
	return transpose(pred)
end

# lors

"""
	lor_maxq(hitdf::DataFrame)

Compute lors using the SiPM of max charge (maxq algo) to divide the event
in two hemispheres and estimating the lor vertex from barycenter
"""
function lor_maxq(hitdf::DataFrame)
	hq2df, hq1df = sipmsel(hitdf)   # select using maxq
	b1 = baricenter(hq1df)          # baricenters
	b2 = baricenter(hq2df)
	return b1, b2, hq1df, hq2df
end


"""
	lor_kmeans(hitdf::DataFrame)

Compute lors using kmeans clustering.
Returns two estimations of vertices. One based in pure kmeans, the
other in barycenter.
"""
function lor_kmeans(hitdf::DataFrame)
	Mhits = get_hits_as_matrix(hitdf)  # take the underlying matrix
	kr = kmeans(Mhits, 2)              # apply kmeans
	ka = assignments(kr) # get the assignments of points to clusters
	kc = counts(kr) # get the cluster sizes
	#rk = kr.centers # get the cluster centers

	hq2df, hq1df = ksipmsel(hitdf, ka)   # select using kmeans list
	b1 = baricenter(hq1df)     # baricenters
	b2 = baricenter(hq2df)
	return b1, b2, hq1df, hq2df
end


"""
	radial_correction(b::Hit, r::Float32, rsipm::Float32)

Take the estimated radius of interaction (r), and the radius of the sipm
and return the corrected positions
"""
function radial_correction(b::Hit, r::Float32, rsipm::Float32)
	r2 = r / rsipm   # redial correction
	return r2 * b.x, r2 * b.y, b.z
end


"""
	phistd(hitdf::DataFrame)

Compute the std deviation in phi weighted by charge, e.g:
Sqrt(1/Q Sum_i (phi_i - phi_mean) * qi )
"""
function phistd(hitdf::DataFrame)
	phi = fphi(hitdf)
	return wstd(phi, hitdf.q)
end


"""
	xyzstd(hitdf::DataFrame)

Compute the std deviation in x weighted by charge, e.g:
Sqrt(1/Q Sum_i (phi_i - phi_mean) * qi )
"""
function xyzstd(hitdf::DataFrame, column::String="x")
	x = hitdf[!, column]
	return wstd(x, hitdf.q)
end


#q correction

"""
	qcor(df::DataFrame, lf::Function, yref::Number)

Linear correction to q
"""
function qcor(df::DataFrame, lf::Function, yref::Number)
    yold = df.y_mean
    ypre = lf.(df.x_mean)
    return (yold .- ypre) .+ yref
end


function qcor2!(df::DataFrame, lf::Function,
	            xc="r", yc="q1", new="qc", yref=2000.0)
    yold = df[!, yc]
    ypre = lf.(df[!, xc])
    df[!, new] = (yold .- ypre) .+ yref
end
