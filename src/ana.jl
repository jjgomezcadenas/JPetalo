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

# """
# 	recohits(event            ::Integer,
# 				   total_charge::DataFrame,
# 				   sensor_xyz  ::DataFrame,
# 				   waveform    ::DataFrame,
# 				   ecut        ::Float32,
# 				   pde         ::Float32,
# 				   max_pes     ::Float32,
# 				   sigma_tof   ::Float32)
#
# For each event, take the three tables (total_charge, sensor_xyz, waveform)
# and produce a xyzqt data frame, where each SiPM has information of position
# (x,y,z), total charge in the SiPM and time (first photon).
#
# """
# function recohits(event        ::Integer,
# 				   total_charge::DataFrame,
# 				   sensor_xyz  ::DataFrame,
# 				   waveform    ::DataFrame,
# 				   ecut        ::Float32,
# 				   pde         ::Float32,
# 				   max_pes     ::Integer,
# 				   sigma_tof   ::Float32)
#
#
# 	function sync!(wfm, qdf)
# 	  j=1
# 	  for i in 1:nrow(wfm)
# 		  if wfm[i, "sensor_id"] ∈ qdf.sensor_id
# 			  nomatch=true
# 			  while nomatch && j <= nrow(qdf)
# 				  if wfm[i, "sensor_id"] == qdf[j, "sensor_id"]
# 					  @debug "wfm sensor_id " wfm[i, "sensor_id"]
# 					  @debug "wfm q_sum " wfm[i, "q_sum"]
# 					  @debug "qdf sensor_id " qdf[j, "sensor_id"]
# 					  @debug "qdf charge " qdf[j, "charge"]
#
# 					  wfm[i, "q_sum"] = qdf[j, "charge"] * pde
# 					  j+=1
# 					  nomatch=false
# 				  else
# 					  j+=1
# 				  end
# 			  end
# 		  end
# 	  end
# 	end
#
#
#   	# select the waveform of this event
# 	wfm = select_by_column_value(waveform, "event_id", event)
# 	if nrow(wfm) == 0
# 		return nothing
# 	end
#
# 	# add a column with probability of surviving pdf cut (pass if prob < pde)
# 	wfm[!, "prob"] = rand(Float32, nrow(wfm))
#
# 	# add a column of charge (each photon arriving to the SiPm in the waveform has charge of one)
# 	wfm[!, "q"]    = Float32.(ones(nrow(wfm)))
#
# 	@debug "event =$event, waveform dataframe: size =$(size(wfm))"
# 	@debug first(wfm, 5)
#
# 	# SiPM pass PDE cut if prob < PDE
# 	wfmc = select_by_column_value_lt(wfm, "prob", pde)
# 	if nrow(wfmc) == 0
# 		return nothing
# 	end
#
# 	# add a column of gaussian random numbers representing the smearing of electronics and sensor
# 	d    = Normal(0.0, sigma_tof)
# 	wfmc[!,"dt"] = Float32.(rand(d, nrow(wfmc)))
#
# 	# add column of smeared times to true times
# 	wfmt = transform(wfmc, [:time, :dt] => (+) => :mtime)
#
# 	@debug "waveform after prob cut and smeared time: size =$(size(wfmt))"
# 	@debug first(wfmt, 5)
#
# 	# Find the resolution (average difference beteen t1 and tr1)
# 	trdt = mean(wfmt.time - wfmt.mtime)
#
# 	# group by SiPMs and take minimum time
# 	wtmin = combine(groupby(wfmt, :sensor_id), :time => minimum)
#
# 	@debug "waveform after grouping SiPMs and min time: size = $(size(wtmin))"
# 	@debug first(wtmin, 5)
#
# 	# group by SiPMs and take minimum time
# 	wrtmin = combine(groupby(wfmt, :sensor_id), :mtime => minimum)
# 	@debug "waveform after grouping SiPMs and min reco time: size = $(size(wrtmin))"
# 	@debug first(wrtmin, 5)
#
# 	trmdt = mean(wtmin.time_minimum - wrtmin.mtime_minimum)
#
# 	# group by SiPMs and compute the sum of the charge in the SiPMs
# 	wtmq = combine(groupby(wfmt, :sensor_id), :q => sum)
#
# 	@debug "waveform after grouping SiPMs and sum charge: size = $(size(wtmq))"
# 	@debug first(wtmq, 5)
#
# 	#select the dataframe of total charge for this event
# 	qdf = select_by_column_value(total_charge, "event_id", event)
# 	if nrow(qdf) == 0
# 		return nothing
# 	end
#
# 	@debug "total charge DF: size = $(size(qdf))"
# 	@debug first(qdf, 5)
#
# 	# find SiPMs which have more than max_pes (max number stored in waveform)
# 	qdf10 = select_by_column_value_gt(qdf, "charge", max_pes)
#
# 	@debug "select sipms such that q > 10 pes: size = $(size(qdf10))"
# 	@debug first(qdf10, 5)
#
# 	#sync the DataFrames: copy to wtmq SiPMs with charge > max_pes
# 	if nrow(qdf10) > 0
# 		sync!(wtmq, qdf10)
# 	end
#
# 	@debug "data frame after sync: size = $(size(wtmq))"
# 	@debug first(wtmq, 5)
#
# 	#construct qt dataframe
# 	wfmx  = DataFrame(sensor_id=wtmq.sensor_id,
# 					   tmin=wtmin.time_minimum,
# 					   trmin=wrtmin.mtime_minimum,
# 					   q=wtmq.q_sum)
#
#   	@debug "DF with time and q: size = $(size(wfmx))"
# 	@debug first(wfmx, 5)
#
# 	# cut on total charge (> ecut)
# 	qdft   = wfmx[wfmx.q .>ecut,:]
#
# 	#select the id's of the sipms with charge above threshold
# 	sids = qdft[!,:sensor_id]
#
# 	#compute positions of the SiPMs
# 	pos   = sipm_pos.((sensor_xyz,),sids)
# 	x = [p[1] for p in pos]
# 	y = [p[2] for p in pos]
# 	z = [p[3] for p in pos]
#
# 	# Construct hit data frame
# 	xyzqt   = DataFrame(x=x,y=y,z=z,
# 					  tmin=qdft.tmin,
# 					  trmin=qdft.trmin,
# 					  q=qdft.q)
#
# 	return xyzqt
# end
