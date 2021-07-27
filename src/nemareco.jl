using DataFrames
using Glob
using Maybe


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
				   max_pes     ::Integer,
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
					  @debug "qdf sensor_id " qdf[j, "sensor_id"]
					  @debug "qdf charge " qdf[j, "charge"]

					  wfm[i, "q_sum"] = qdf[j, "charge"] * pde
					  j+=1
					  nomatch=false
				  else
					  j+=1
				  end
			  end
		  end
	  end
	end


  	# select the waveform of this event
	wfm = select_by_column_value(waveform, "event_id", event)
	if nrow(wfm) == 0
		return nothing
	end

	# add a column with probability of surviving pdf cut (pass if prob < pde)
	wfm[!, "prob"] = rand(Float32, nrow(wfm))

	# add a column of charge (each photon arriving to the SiPm in the waveform has charge of one)
	wfm[!, "q"]    = Float32.(ones(nrow(wfm)))

	@debug "event =$event, waveform dataframe: size =$(size(wfm))"
	@debug first(wfm, 5)

	# SiPM pass PDE cut if prob < PDE
	wfmc = select_by_column_value_lt(wfm, "prob", pde)
	if nrow(wfmc) == 0
		return nothing
	end

	# add a column of gaussian random numbers representing the smearing of electronics and sensor
	d    = Normal(0.0, sigma_tof)
	wfmc[!,"dt"] = Float32.(rand(d, nrow(wfmc)))

	# add column of smeared times to true times
	wfmt = transform(wfmc, [:time, :dt] => (+) => :mtime)

	@debug "waveform after prob cut and smeared time: size =$(size(wfmt))"
	@debug first(wfmt, 5)

	# Find the resolution (average difference beteen t1 and tr1)
	trdt = mean(wfmt.time - wfmt.mtime)

	# group by SiPMs and take minimum time
	wtmin = combine(groupby(wfmt, :sensor_id), :time => minimum)

	@debug "waveform after grouping SiPMs and min time: size = $(size(wtmin))"
	@debug first(wtmin, 5)

	# group by SiPMs and take minimum time
	wrtmin = combine(groupby(wfmt, :sensor_id), :mtime => minimum)
	@debug "waveform after grouping SiPMs and min reco time: size = $(size(wrtmin))"
	@debug first(wrtmin, 5)

	trmdt = mean(wtmin.time_minimum - wrtmin.mtime_minimum)

	# group by SiPMs and compute the sum of the charge in the SiPMs
	wtmq = combine(groupby(wfmt, :sensor_id), :q => sum)

	@debug "waveform after grouping SiPMs and sum charge: size = $(size(wtmq))"
	@debug first(wtmq, 5)

	#select the dataframe of total charge for this event
	qdf = select_by_column_value(total_charge, "event_id", event)
	if nrow(qdf) == 0
		return nothing
	end

	@debug "total charge DF: size = $(size(qdf))"
	@debug first(qdf, 5)

	# find SiPMs which have more than max_pes (max number stored in waveform)
	qdf10 = select_by_column_value_gt(qdf, "charge", max_pes)

	@debug "select sipms such that q > 10 pes: size = $(size(qdf10))"
	@debug first(qdf10, 5)

	#sync the DataFrames: copy to wtmq SiPMs with charge > max_pes
	if nrow(qdf10) > 0
		sync!(wtmq, qdf10)
	end

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
						rsipm       ::Float32,
						detconf     ::String,
						df1         ::DataFrame,
						df2         ::DataFrame,
						primaries   ::DataFrame,
						total_charge::DataFrame,
						sensor_xyz  ::DataFrame,
						waveform    ::DataFrame,
						lor_algo    ::Function,
						n3d         ::Dict)

Extracts the information needed for nema studies, including lors
"""
function nema_analysis!(event       ::Integer,
	                    rsipm       ::Float32,
						dc          ::DetConf,
						df1         ::DataFrame,
						df2         ::DataFrame,
						primaries   ::DataFrame,
						total_charge::DataFrame,
						sensor_xyz  ::DataFrame,
						waveform    ::DataFrame,
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
	prim = select_by_column_value(primaries, "event_id", event)
	@debug "Primaries in event:" prim
	#hit dataframe
	hitdf = recohits(event, total_charge, sensor_xyz, waveform,
					  dc.ecut, dc.pde, dc.max_pes, dc.sigma_tof)

	if hitdf == nothing
		@warn "Warning, hidtf evaluates to nothing for event = $event"
		return
	end

	if nrow(hitdf) < 2
		@warn "Warning, hidtf is <2 for event = $event"
		return
	end

	@info " hit dataframe: size = $size(hitdf)"
	@debug first(hitdf, 5)

	# reconstruct (x,y) : barycenter
	#b1, b2, hq1df, hq2df = lor_kmeans(hitdf)
	#b1, b2, hq1df, hq2df = lor_maxq(hitdf)
	b1, b2, hq1df, hq2df = lor_algo(hitdf)

	@info " barycenters" b1 b2

	# total charge
	q1 = sum(hq1df.q)
	q2 = sum(hq2df.q)

	@info " total charge: q1 = $q1, q2 = $q2"
	if q1 < dc.qmin || q1 > dc.qmax
		@warn "Warning, q1 is $q1 for event $event"
		return false
	end
	if q2 < dc.qmin || q2 > dc.qmax
		@warn "Warning, q2 is $q2 for event $event"
		return false
	end

	# Compute phistd and zstd1
	#phistd1 = phistd(hq1df)
	#zstd1   = xyzstd(hq1df,"z")
	#phistd2 = phistd(hq2df)
	#zstd2   = xyzstd(hq2df,"z")
	#@info " phistd1 = $phistd1, zstd1 = $zstd1"
	#@info " phistd2 = $phistd2, zstd2 = $zstd2"

	# find true position (and correlate with barycenter)
	xt1, xt2 = true_xyz(b1, b2, df1, df2)
	# find r1 and r2 (from True info)
	r1 = rxy(xt1[1], xt1[2])
	r2 = rxy(xt2[1], xt2[2])
	@info " True position in hemisphere 1" xt1
	@info " True position in hemisphere 1" xt2

	# find r1 and r2 from charge
	r1q = dc.rfq(q1)
	r2q = dc.rfq(q2)

	# find r1 and r2 from zstd
	#r1z = dc.rfz(zstd1)
	#r2z = dc.rfz(zstd2)
	#@info " True     : r1 = $r1, r2 = $r2"
	#@info " From r   : r1q = $r1q, r2q = $r2q"
	#@info " From zstd: r1z = $r1z, r2z = $r2z"

	# New (x,y) positions estimated from r1, r2
	x1, y1, z1  = radial_correction(b1, r1, rsipm)
	x2, y2, z2  = radial_correction(b2, r2, rsipm)

	# New (x,y) positions estimated from r1q, r2q
	xr1, yr1, zr1  = radial_correction(b1, r1q, rsipm)
	xr2, yr2, zr2  = radial_correction(b2, r2q, rsipm)

	# New (x,y) positions estimated from r1z, r2z
	#xR1, yR1, zR1  = radial_correction(b1, r1z, rsipm)
	#xR2, yR2, zR2  = radial_correction(b2, r2z, rsipm)

	@info " New (x,y,z) positions estimated from r1, r2 & r1q, r2q"
	@info " from r1:  x1 = $x1, y1=$y1, z1=$z1"
	@info " from r2:  x2 = $x2, y1=$y2, z1=$z2"
	@info " from rq: xr1 = $xr1, yr1=$yr1, zr1=$zr1"
	@info " from rq: xr2 = $xr2, y1=$yr2, z1=$zr2"
	#@info " from rz: xR1 = $xR1, yR1=$yR1, zR1=$zR1"
	#@info " from rz: xR2 = $xR2, yR2=$yR2, zR2=$zR2"


	# Find the sipm with the fastest time
	t1 = minimum(hq1df.tmin)
	t2 = minimum(hq2df.tmin)
	tr1 = minimum(hq1df.trmin)
	tr2 = minimum(hq2df.trmin)

	ntof1 = min(dc.ntof, nrow(hq1df))
	ntof2 = min(dc.ntof, nrow(hq2df))

	# sort reco times in ascending order
	t1s = sort(hq1df.trmin)
	t2s = sort(hq2df.trmin)
	# take average
	ta1 = mean(t1s[1:ntof1])
	ta2 = mean(t2s[1:ntof2])

	@info " time"
	@info " true:  t1 = $t1, t2=$t2"
	@info " smeared:  tr1 = $tr1, tr2=$tr2"
	@info " averaged:  ta1 = $ta1, ta2=$ta2"

	ht1 = select_by_column_value(hq1df, "tmin", t1)
	ht2 = select_by_column_value(hq2df, "tmin", t2)
	htr1 = select_by_column_value(hq1df, "trmin", tr1)
	htr2 = select_by_column_value(hq2df, "trmin", tr2)
	# position of the SiPM with fastet signal.
	# xct1 = ht1.x[1]
	# yct1 = ht1.y[1]
	# zct1 = ht1.z[1]
	# xct2 = ht2.x[1]
	# yct2 = ht2.y[1]
	# zct2 = ht2.z[1]


	# store data
	#source data

	push!(n3d["xs"],prim.x[1])
	push!(n3d["ys"],prim.y[1])
	push!(n3d["zs"],prim.z[1])
	push!(n3d["ux"],prim.vx[1])
	push!(n3d["uy"],prim.vy[1])
	push!(n3d["uz"],prim.vz[1])
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
	push!(n3d["t1"],t1)
	push!(n3d["x2"],x2)
	push!(n3d["y2"],y2)
	push!(n3d["z2"],z2)
	push!(n3d["t2"],t2)
	push!(n3d["xr1"],xr1)
	push!(n3d["yr1"],yr1)
	push!(n3d["zr1"],zr1)
	push!(n3d["tr1"],tr1)
	#push!(n3d["xR1"],xR1)
	#push!(n3d["yR1"],yR1)
	#push!(n3d["zR1"],zR1)
	push!(n3d["xr2"],xr2)
	push!(n3d["yr2"],yr2)
	push!(n3d["zr2"],zr2)
	push!(n3d["tr2"],tr2)
	#push!(n3d["xR2"],xR2)
	#push!(n3d["yR2"],yR2)
	#push!(n3d["zR2"],zR2)
	push!(n3d["ta1"],ta1)
	push!(n3d["ta2"],ta2)
	push!(n3d["xb1"],ht1.x[1])
	push!(n3d["yb1"],ht1.y[1])
	push!(n3d["zb1"],ht1.z[1])
	push!(n3d["xb2"],ht2.x[1])
	push!(n3d["yb2"],ht2.y[1])
	push!(n3d["zb2"],ht2.z[1])

	# single hemisphere data (control)
	push!(n3d["nsipm1"],nrow(hq1df))
	push!(n3d["q1"], sum(hq1df.q))
	push!(n3d["r1"],r1)
	push!(n3d["r1q"],r1q)
	#push!(n3d["r1z"],r1z)
	#push!(n3d["phistd1"],phistd1)
	#push!(n3d["zstd1"],zstd1)

	push!(n3d["nsipm2"],nrow(hq2df))
	push!(n3d["q2"], sum(hq2df.q))
	push!(n3d["r2"],r2)
	push!(n3d["r2q"],r2q)
	#push!(n3d["r2z"],r2z)
	#push!(n3d["phistd2"],phistd2)
	#push!(n3d["zstd2"],zstd2)

end

"""
	nemareco(files         ::Vector{String},
				  file_i   ::Integer=1,
				  file_l   ::Integer=1,
				  pde      ::Float32=0.3,
				  max_pes  ::Integer=10,
				  sigma_tof::Float32=0.085,
				  ecut     ::Float32=2.0,
				  qmin     ::Float32=1.400,
				  qmax     ::Float32=3000.0,
				  cq1      ::Float32=297.9,
				  cq2      ::Float32=0.0031,
				  ntof     ::Integer = 5,
				  phot     ::Bool=true,
				  lor_algo ::Function=lor_maxq)

Main driver for nema studies

"""
function nemareco(files    ::Vector{String},
				  dconf    ::DetConf,
	              file_i   ::Integer=1,
				  file_l   ::Integer=1,
				  phot     ::Bool=true,
				  lor_algo ::Function=lor_maxq)

	# define data dictionary

	n3d = Dict("nsipm1"=>[0],"nsipm2"=>[0],
			   "q1" =>[0.0f0],   "q2" =>[0.0f0],
	           "r1"  =>[0.0f0],  "r2"  =>[0.0f0],
			   "r1q"  =>[0.0f0], "r2q"  =>[0.0f0],
			   #"phistd1"=>[0.0f0],  "zstd1"=>[0.0f0],
			    #"r1z"  =>[0.0f0],
			   #"phistd2"=>[0.0f0],  "zstd2"=>[0.0f0],
		   	   #"r2z"  =>[0.0f0],
			   "xs"=>[0.0f0], "ys"=>[0.0f0], "zs"=>[0.0f0],
		       "ux"=>[0.0f0], "uy"=>[0.0f0], "uz"=>[0.0f0],
	           "xt1"=>[0.0f0], "yt1"=>[0.0f0], "zt1"=>[0.0f0],
		       "xt2"=>[0.0f0], "yt2"=>[0.0f0], "zt2"=>[0.0f0],
               "x1"=>[0.0f0],   "y1"=>[0.0f0], "z1"=>[0.0f0],"t1"=>[0.0f0],
               "x2"=>[0.0f0],   "y2"=>[0.0f0], "z2"=>[0.0f0], "t2"=>[0.0f0],
			   "xr1"=>[0.0f0], "yr1"=>[0.0f0], "zr1"=>[0.0f0], "tr1"=>[0.0f0],
			   #"xR1"=>[0.0f0], "yR1"=>[0.0f0], "zR1"=>[0.0f0],
               "xr2"=>[0.0f0], "yr2"=>[0.0f0], "zr2"=>[0.0f0], "tr2"=>[0.0f0],
			   #"xR2"=>[0.0f0], "yR2"=>[0.0f0], "zR2"=>[0.0f0],
			   "xb1"=>[0.0f0], "yb1"=>[0.0f0], "zb1"=>[0.0f0],
			   "xb2"=>[0.0f0], "yb2"=>[0.0f0], "zb2"=>[0.0f0],
			   "ta1"=>[0.0f0], "ta2"=>[0.0f0])

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
						nema_analysis!(event, rsipm, dconf,
						               df1, df2, pdf.primaries,
									   pdf.total_charge, pdf.sensor_xyz, pdf.waveform,
									   lor_algo, n3d)

        			end
				else
					nema_analysis!(event, rsipm, dconf,
								   df1, df2, pdf.primaries,
								   pdf.total_charge, pdf.sensor_xyz, pdf.waveform,
								   lor_algo, n3d)
				end
			end
    	end
	end
	n3df = DataFrame(n3d)
	return n3df[2:end,:]
end
