using DataFrames
using Glob


function grfq(q)
    return Float32(297.9 + 0.031 * q)
end
rfq(q::Float32) = grfq(q)
rfq(q::Float64) = grfq(q)


function nema3a(files::Vector{String}, file_i::Integer=1, file_l::Integer=1,
				ecut::Float32=2.0, pde::Float32=0.3,
				qmin::Float32=1.400, qmax::Float32=3000.0,
				phot=true, warn=false)

	function qlimit(q::Float32, qmin::Float32, qmax::Float32, event::Integer)
		if q < qmin || q > qmax
			if warn
				println("Warning, q is ", q, " for event ", event)
			end
			return false
		end
	end

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

	function photo_analysis!(event::Integer,
		                    df1::DataFrame, df2::DataFrame,
							qdf::DataFrame, sxyzdf::DataFrame)

		#select reco
		hitdf  = reco_hits(event, ecut, pde, qdf, sxyzdf)

		if nrow(hitdf) < 2
			if warn
				println("Warning, hidtf is <2 for event =", event)
			end
			return
		end

		# reconstruct (x,y) : barycenter
		#println(event, nrow(hitdf))
		b1, b2, hq1df, hq2df = lor_kmeans(hitdf)

		# total charge
		q1 = sum(hq1df.q)
		q2 = sum(hq2df.q)

		if qlimit(q1, qmin, qmax, event) == false
			return
		end
		if qlimit(q2, qmin, qmax, event) == false
			return
		end

		# Compute phistd and zstd1
		phistd1 = phistd(hq1df)
		zstd1 = xyzstd(hq1df,"z")

		# find true position (and correlate with barycenter)
		xt1, xt2 = true_xyz(b1, b2, df1, df2)
		# find r1 and r2 (from True info)
		r1 = rxy(xt1[1], xt1[2])
        r2 = rxy(xt2[1], xt2[2])

		# find r1 and r2 from charge
		r1q = rfq(q1)
		r2q = rfq(q2)

		# New (x,y) positions estimated from r1, r2
		x1, y1, z1  = radial_correction(b1, r1, rsipm)
		x2, y2, z2  = radial_correction(b2, r2, rsipm)

		# New (x,y) positions estimated from r1q, r2q
		xr1, yr1, zr1  = radial_correction(b1, r1q, rsipm)
		xr2, yr2, zr2  = radial_correction(b2, r2q, rsipm)

		# store data
		push!(n3d["nsipm"],nrow(hq1df))
		xyzt!(xt1, xt2)
		#xyzb!(b1,b2)
		xyzr!(x1,y1,z1,x2,y2,z2)
	  	xyzrq!(xr1,yr1,zr1,xr2,yr2,zr2)
		xyzq!(hq1df)
		phiq!(phistd1, zstd1)
		r12!(r1, r1q)
	end

	function xyzt!(xt1::Vector{Float32}, xt2::Vector{Float32})
		push!(n3d["xt1"],xt1[1])
		push!(n3d["yt1"],xt1[2])
		push!(n3d["zt1"],xt1[3])
		push!(n3d["xt2"],xt2[1])
		push!(n3d["yt2"],xt2[2])
		push!(n3d["zt2"],xt2[3])
	end

	function xyzb!(b1::Hit, b2::Hit)
		push!(n3d["xb1"],b1.x)
		push!(n3d["yb1"],b1.y)
		push!(n3d["zb1"],b1.z)
		push!(n3d["xb2"],b2.x)
		push!(n3d["yb2"],b2.y)
		push!(n3d["zb2"],b2.z)
	end

	function xyzr!(x1::Float32, y1::Float32, z1::Float32,
		           x2::Float32, y2::Float32, z2::Float32)

		push!(n3d["x1"],x1)
		push!(n3d["y1"],y1)
		push!(n3d["z1"],z1)
		push!(n3d["x2"],x2)
		push!(n3d["y2"],y2)
		push!(n3d["z2"],z2)
	end

	function xyzrq!(x1::Float32, y1::Float32, z1::Float32,
		           x2::Float32, y2::Float32, z2::Float32)
		push!(n3d["xr1"],x1)
		push!(n3d["yr1"],y1)
		push!(n3d["zr1"],z1)
		push!(n3d["xr2"],x2)
		push!(n3d["yr2"],y2)
		push!(n3d["zr2"],z2)
end

	function xyzq!(hdf::DataFrame)
		#push!(n3d["qmx"], maximum(hdf.q))
		push!(n3d["q1"], sum(hdf.q))
		#push!(n3d["qav"], mean(hdf.q))
		#push!(n3d["dz"], maximum(abs.(hdf.z)) - minimum(abs.(hdf.z)))
	end

	function phiq!(phistd::Float32, zstd::Float32)
		push!(n3d["phistd"],phistd)
		push!(n3d["zstd"],zstd)
	end

	function r12!(r::Float32, rq::Float32)
		push!(n3d["r"],r)
		push!(n3d["rq"],rq)
	end

	# define data dictionary
	n3d = Dict("nsipm"=>[0],
	           "r"  =>[0.0], "phistd"=>[0.0],  "zstd"=>[0.0],
			   "rq"  =>[0.0], "q1" =>[0.0],
			   #"qmx"=>[0.0],  "qav"=>[0.0],
	           "xt1"=>[0.0], "yt1"=>[0.0], "zt1"=>[0.0],
		       "xt2"=>[0.0], "yt2"=>[0.0], "zt2"=>[0.0],
               "x1"=>[0.0],   "y1"=>[0.0], "z1"=>[0.0],
               "x2"=>[0.0],   "y2"=>[0.0], "z2"=>[0.0],
			   "xr1"=>[0.0], "yr1"=>[0.0], "zr1"=>[0.0],
               "xr2"=>[0.0], "yr2"=>[0.0], "zr2"=>[0.0])
			   #"xb1"=>[0.0], "yb1"=>[0.0], "zb1"=>[0.0],
			   #"xb2"=>[0.0], "yb2"=>[0.0], "zb2"=>[0.0])

	# read one file to compute the radius of sipm
	pdf = read_abc(files[1])
	rsipm = rxy(pdf.sensor_xyz.x[2], pdf.sensor_xyz.y[2])
	#sqrt(pdf.sensor_xyz.x[2]^2+pdf.sensor_xyz.y[2]^2)

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
						photo_analysis!(event, df1, df2,
						               pdf.total_charge, pdf.sensor_xyz)
        			end
				else
						photo_analysis!(event, df1, df2,
									   pdf.total_charge, pdf.sensor_xyz)
				end
			end
    	end
	end
	n3df = DataFrame(n3d)
	return n3df[2:end,:]
end
