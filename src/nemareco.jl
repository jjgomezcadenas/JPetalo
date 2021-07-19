using DataFrames
using Glob

"""
	rfq(q)
radius-from-q function. The parameters of the straight line have been
obtained from a fit to photoelectric data.
"""
function grfq(q)
    return Float32(326.9 + 0.0226 * q)
end
rfq(q::Float32) = grfq(q)
rfq(q::Float64) = grfq(q)


"""
	recohits(event            ::Integer,
				   total_charge::DataFrame,
				   sensor_xyz  ::DataFrame,
				   waveform    ::DataFrame,
				   ecut        ::Float32,
				   pde         ::Float32,
				   max_pes     ::Float32,
				   sigma_tof   ::Float32)

For each event, take the three tables (total_charge, sensor_xyz, waveform)
and produce a xyzqt data frame, where each SiPM has information of position
(x,y,z), total charge in the SiPM and time (first photon).

"""
function recohits(event        ::Integer,
				   total_charge::DataFrame,
				   sensor_xyz  ::DataFrame,
				   waveform    ::DataFrame,
				   ecut        ::Float32,
				   pde         ::Float32,
				   max_pes     ::Float32,
				   sigma_tof   ::Float32)


	function sync!(wfm, qdf)
	  j=1
	  for i in 1:nrow(wfm)
		  if wfm[i, "sensor_id"] ∈ qdf.sensor_id
			  nomatch=true
			  while nomatch && j <= nrow(qdf)
				  if wfm[i, "sensor_id"] == qdf[j, "sensor_id"]
					  @debug "wfm sensor_id " wfm[i, "sensor_id"]
					  @debug "wfm q_sum " wfm[i, "q_sum"]
					  @debug "qdf sensor_id " qdf[i, "sensor_id"]
					  @debug "qdf charge " qdf[i, "charge"]

					  wfm[i, "q_sum"] = qdf[j, "charge"] * pde
					  j+=1
					  nomatch=false
				  else
					  j+=2
				  end
			  end
		  end
	  end
	end


  	# select the waveform of this event
	wfm = select_by_column_value(waveform, "event_id", event)

	# add a column with probability of surviving pdf cut (pass if prob < pde)
	wfm[!, "prob"] = rand(Float32, nrow(wfm))

	# add a column of charge (each photon arriving to the SiPm in the waveform has charge of one)
	wfm[!, "q"]    = Float32.(ones(nrow(wfm)))

	@debug "event =$event, waveform dataframe: size =$(size(wfm))"
	@debug first(wfm, 5)

	# SiPM pass PDE cut if prob < PDE
	wfmc = select_by_column_value_lt(wfm, "prob", pde)

	# add a column of gaussian random numbers representing the smearing of electronics and sensor
	d    = Normal(0.0, sigma_tof)
	wfmc[!,"dt"] = Float32.(rand(d, nrow(wfmc)))

	# add column of smeared times to true times
	wfmt = transform(wfmc, [:time, :dt] => (+) => :mtime)

	@debug "waveform after prob cut and smeared time: size =$(size(wfmt))"
	@debug first(wfmt, 5)

	# group by SiPMs and take minimum time
	wtmin = combine(groupby(wfmt, :sensor_id), :time => minimum)

	@debug "waveform after grouping SiPMs and min time: size = $(size(wtmin))"
	@debug first(wtmin, 5)

	# group by SiPMs and take minimum time
	wrtmin = combine(groupby(wfmt, :sensor_id), :mtime => minimum)
	@debug "waveform after grouping SiPMs and min reco time: size = $(size(wrtmin))"
	@debug first(wrtmin, 5)

	# group by SiPMs and compute the sum of the charge in the SiPMs
	wtmq = combine(groupby(wfmt, :sensor_id), :q => sum)

	@debug "waveform after grouping SiPMs and sum charge: size = $(size(wtmq))"
	@debug first(wtmq, 5)

	#select the dataframe of total charge for this event
	qdf = select_by_column_value(total_charge, "event_id", event)

	@debug "total charge DF: size = $(size(qdf))"
	@debug first(qdf, 5)

	# find SiPMs which have more than max_pes (max number stored in waveform)
	qdf10 = select_by_column_value_gt(qdf, "charge", max_pes)

	@debug "select sipms such that q > 10 pes: size = $(size(qdf10))"
	@debug first(qdf10, 5)

	#sync the DataFrames: copy to wtmq SiPMs with charge > max_pes
	sync!(wtmq, qdf10)

	@debug "data frame after sync: size = $(size(wtmq))"
	@debug first(wtmq, 5)

	#construct qt dataframe
	  wfmx  = DataFrame(sensor_id=wtmq.sensor_id,
					   tmin=wtmin.time_minimum,
					   trmin=wrtmin.mtime_minimum,
					   q=wtmq.q_sum)

  	@debug "DF with time and q: size = $(size(wfmx))"
	@debug first(wfmx, 5)

	# cut on total charge (> ecut)
	qdft   = wfmx[wfmx.q .>ecut,:]

	#select the id's of the sipms with charge above threshold
	sids = qdft[!,:sensor_id]

	#compute positions of the SiPMs
	pos   = sipm_pos.((sensor_xyz,),sids)
	x = [p[1] for p in pos]
	y = [p[2] for p in pos]
	z = [p[3] for p in pos]

	# Construct hit data frame
	xyzqt   = DataFrame(x=x,y=y,z=z,
					  tmin=qdft.tmin,
					  trmin=qdft.trmin,
					  q=qdft.q)

	return xyzqt
end

"""
	nema_analysis!(event       ::Integer,
							 df1         ::DataFrame,
							 df2         ::DataFrame,
							 total_charge::DataFrame,
							 sensor_xyz  ::DataFrame,
							 waveform    ::DataFrame,
							 ecut        ::Float32,
							 pde         ::Float32,
							 max_pes     ::Float32,
							 sigma_tof   ::Float32,
							 lor_algo    ::Function)

Extracts the information needed for nema studies, including lors
"""
function nema_analysis!(event       ::Integer,
						 df1         ::DataFrame,
						 df2         ::DataFrame,
						 primaries   ::DataFrame,
						 total_charge::DataFrame,
						 sensor_xyz  ::DataFrame,
						 waveform    ::DataFrame,
						 ecut        ::Float32,
						 pde         ::Float32,
						 max_pes     ::Integer,
						 sigma_tof   ::Float32,
						 lor_algo    ::Function,
						 n3d         ::Dict)

	# true_xyz function
	function true_xyz(b1::Hit, b2::Hit, df1::DataFrame, df2::DataFrame)
		# distance between baricenter and true
		d1 = dxyz([b1.x, b1.y, b1.z], [df1.x[1], df1.y[1], df1.z[1]])
		d2 = dxyz([b1.x, b1.y, b1.z], [df2.x[1], df2.y[1], df2.z[1]])

		if d2 < d1
			 xt2 = [df1.x[1], df1.y[1], df1.z[1]]
			 xt1 = [df2.x[1], df2.y[1], df2.z[1]]
		 else
			 xt1 = [df1.x[1], df1.y[1], df1.z[1]]
			 xt2 = [df2.x[1], df2.y[1], df2.z[1]]
		 end
		 return xt1, xt2
	end

	# primaries
	pdf = select_by_column_value(primaries, "event_id", event)
	#hit dataframe
	hitdf  = recohits(event, total_charge, sensor_xyz, waveform,
					  ecut, pde, max_pes, sigma_tof)

	if nrow(hitdf) < 2
		@warn "Warning, hidtf is <2 for event = $event"
		return
	end

	@info " hit dataframe: size = $size(hitdf)"
	@info first(hitdf, 5)

	# reconstruct (x,y) : barycenter
	#b1, b2, hq1df, hq2df = lor_kmeans(hitdf)
	#b1, b2, hq1df, hq2df = lor_maxq(hitdf)
	b1, b2, hq1df, hq2df = lor_algo(hitdf)

	@info " barycenters" b1 b2

	# total charge
	q1 = sum(hq1df.q)
	q2 = sum(hq2df.q)

	@info " total charge: q1 = $q1, q2 = $q2"
	if q1 < qmin || q1 > qmax
		@warn "Warning, q1 is $q1 for event $event$"
		return false
	end
	if q2 < qmin || q2 > qmax
		@warn "Warning, q2 is $q2 for event $event$"
		return false
	end

	# Compute phistd and zstd1
	phistd1 = phistd(hq1df)
	zstd1 = xyzstd(hq1df,"z")
	@info " phistd1 = $phistd1, zstd1 = $zstd1"

	# find true position (and correlate with barycenter)
	xt1, xt2 = true_xyz(b1, b2, df1, df2)
	# find r1 and r2 (from True info)
	r1 = rxy(xt1[1], xt1[2])
	r2 = rxy(xt2[1], xt2[2])
	@info " True position in hemisphere 1" xt1
	@info " True position in hemisphere 1" xt2

	# find r1 and r2 from charge
	r1q = rfq(q1)
	r2q = rfq(q2)
	@info " True: r1 = $r1, r2 = $r2"
	@info " True: r1q = $r1q, r2q = $r2q"

	# New (x,y) positions estimated from r1, r2
	x1, y1, z1  = radial_correction(b1, r1, rsipm)
	x2, y2, z2  = radial_correction(b2, r2, rsipm)

	# New (x,y) positions estimated from r1q, r2q
	xr1, yr1, zr1  = radial_correction(b1, r1q, rsipm)
	xr2, yr2, zr2  = radial_correction(b2, r2q, rsipm)

	@info " New (x,y,z) positions estimated from r1, r2 & r1q, r2q"
	@info " from r1:  x1 = $x1, y1=$y1, z1=$z1"
	@info " from r2:  x2 = $x2, y1=$y2, z1=$z2"
	@info " from rb1: xr1 = $xr1, yr1=$yr1, zr1=$zr1"
	@info " from rb2: xr2 = $xr2, y1=$yr2, z1=$zr2"

	t1 = hq1df.tmin
	t2 = hq2df.tmin

	# store data
	#source data
	push!(n3d["xs"],pdf.x)
	push!(n3d["ys"],pdf.y)
	push!(n3d["zs"],pdf.z)
	push!(n3d["ux"],pdf.vx)
	push!(n3d["uy"],pdf.vy)
	push!(n3d["uz"],pdf.vz)
	# double hemisphere data (lors)
	push!(n3d["xt1"],xt1[1])
	push!(n3d["yt1"],xt1[2])
	push!(n3d["zt1"],xt1[3])
	push!(n3d["xt2"],xt2[1])
	push!(n3d["yt2"],xt2[2])
	push!(n3d["zt2"],xt2[3])
	push!(n3d["x1"],x1)
	push!(n3d["y1"],y1)
	push!(n3d["z1"],z1)
	push!(n3d["t1"],hq1df.tmin)
	push!(n3d["x2"],x2)
	push!(n3d["y2"],y2)
	push!(n3d["z2"],z2)
	push!(n3d["t2"],hq2df.tmin)
	push!(n3d["xr1"],x1)
	push!(n3d["yr1"],y1)
	push!(n3d["zr1"],z1)
	push!(n3d["tr1"],hq1df.trmin)
	push!(n3d["xr2"],x2)
	push!(n3d["yr2"],y2)
	push!(n3d["zr2"],z2)
	push!(n3d["tr2"],hq2df.trmin)

	# single hemisphere data (control)
	push!(n3d["nsipm"],nrow(hq1df))
	push!(n3d["q"], sum(hq1df.q))
	push!(n3d["r"],r1)
	push!(n3d["rq"],rq1)
	push!(n3d["phistd"],phistd1)
	push!(n3d["zstd"],zstd1)

end

"""
	nemareco(files    ::Vector{String},
				  file_i   ::Integer=1,
				  file_l   ::Integer=1,
				  pde      ::Float32=0.3,
				  max_pes  ::Integer=10,
				  sigma_tof::Float32=0.085,
				  ecut     ::Float32=2.0,
				  qmin     ::Float32=1.400,
				  qmax     ::Float32=3000.0,
				  prteach  ::Integer=100,
				  phot     ::Bool=true,
				  lor_algo ::Function=lor_algo)

Main driver for nema studies

"""
function nemareco(files    ::Vector{String},
	              file_i   ::Integer=1,
				  file_l   ::Integer=1,
				  pde      ::Float32=0.3,
				  max_pes  ::Integer=10,
                  sigma_tof::Float32=0.085,
				  ecut     ::Float32=2.0,
				  qmin     ::Float32=1.400,
				  qmax     ::Float32=3000.0,
				  prteach  ::Integer=100,
				  phot     ::Bool=true,
				  lor_algo ::Function=lor_maxq)



	# define data dictionary
	n3d = Dict("nsipm"=>[0],
	           "r"  =>[0.0], "phistd"=>[0.0],  "zstd"=>[0.0],
			   "rq"  =>[0.0], "q" =>[0.0],
			   "xs"=>[0.0], "ys"=>[0.0], "zs"=>[0.0],
		       "ux"=>[0.0], "uy"=>[0.0], "ux"=>[0.0],
	           "xt1"=>[0.0], "yt1"=>[0.0], "zt1"=>[0.0],
		       "xt2"=>[0.0], "yt2"=>[0.0], "zt2"=>[0.0],
               "x1"=>[0.0],   "y1"=>[0.0], "z1"=>[0.0],"t1"=>[0.0],
               "x2"=>[0.0],   "y2"=>[0.0], "z2"=>[0.0], "t2"=>[0.0],
			   "xr1"=>[0.0], "yr1"=>[0.0], "zr1"=>[0.0], "tr1"=>[0.0],
               "xr2"=>[0.0], "yr2"=>[0.0], "zr2"=>[0.0], "tr2"=>[0.0])

	# read one file to compute the radius of sipm
	pdf = read_abc(files[1])
	rsipm = rxy(pdf.sensor_xyz.x[2], pdf.sensor_xyz.y[2])

	for file in files[file_i:file_l]               # loop on files
		println("reading file = ", file)
		pdf = read_abc(file)            # read file
		dfs = primary_in_lxe(pdf.vertices)       # primary photons in LXe

		#cevt = -1
		for event in unique(dfs.event_id)       #loop on events
			#  event DF
			vdf = select_by_column_value(dfs, "event_id", event)

			# two primary photons in LXe
			if any(vdf.track_id .== 1) && any(vdf.track_id .== 2)
				df1 = select_by_column_value(vdf, "track_id", 1)
            	df2 = select_by_column_value(vdf, "track_id", 2)
				if phot == true
            		if df1.process_id[1] == 1 && df2.process_id[1] == 1
						nema_analysis!(event,
						               df1, df2, pdf.primaries,
									   pdf.total_charge, pdf.sensor_xyz, pdf.waveform,
									   ecut, pde, max_pes, sigma_tof, lor_algo)

        			end
				else
					nema_analysis!(event,
								   df1, df2, pdf.primaries,
								   pdf.total_charge, pdf.sensor_xyz, pdf.waveform,
								   ecut, pde, max_pes, sigma_tof, lor_algo)
				end
			end
    	end
	end
	n3df = DataFrame(n3d)
	return n3df[2:end,:]
end
