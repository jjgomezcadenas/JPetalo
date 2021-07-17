push!(LOAD_PATH,"../src/")
using DrWatson
@quickactivate "JPetalo"
include("../src/jpetalo.jl")

using DataFrames
using Glob
using ArgParse

"""
	makelor_3p(files::Vector{String},
			   ecut::Float64=4.0, file_i=1, file_l=1)

	Take a vector of files containing PETALO data and compute LORS
	TLR  : A DataFrame of True Lors
	RLHB : A DataFrame of Reco Lors with HALF algorithm and barycenter
	RLKB : A DataFrame of Reco Lors with kmeans algorithm and barycenter
	RLKK : A DataFrame of Reco Lors with kmeans algorithm and kcluster

	Events are pre-selected to be primary photoelectric in LXe
	with two and only two gammas

"""
function makelor_3p(files::Vector{String}, fxy=JPetalo.lor_kmeans,
	                ecut::Float64=4.0, file_i=1, file_l=1)

	# true reconstructed LOR
	TLR = [(t1=0.0,t2=0.0,x1=0.0,y1=0.0,z1=0.0,x2=0.0,y2=0.0,z2=0.0)]
	# LOR reconstructed from barycenter
	RLB = [(t1=0.0,t2=0.0,x1=0.0,y1=0.0,z1=0.0,x2=0.0,y2=0.0,z2=0.0)]
	# LOR corrected using R
	RLR = [(t1=0.0,t2=0.0,x1=0.0,y1=0.0,z1=0.0,x2=0.0,y2=0.0,z2=0.0)]
	dr   = 5 # mm

	pdf = JPetalo.read_abc(files[1])
	rsipm = sqrt(pdf.sensor_xyz.x[2]^2+pdf.sensor_xyz.y[2]^2)
	rmax = rsipm + dr
	rmin = rsipm - dr

	for file in files[file_i:file_l]               # loop on files
		println("reading file = ", file)
		pdf = JPetalo.read_abc(file)            # read file
		dfs = JPetalo.primary_phot_in_lxe(pdf)  # primary photons in LXe

		cevt = 0

		for event in dfs.event_id              #loop on eventa

			#  event DF
			vdf = JPetalo.select_by_column_value(dfs, "event_id", event)
			if nrow(vdf) == 2  && event != cevt    # only 2 gammas
				# select true
				df1 = JPetalo.select_by_column_value(vdf, "track_id", 1)
				df2 = JPetalo.select_by_column_value(vdf, "track_id", 2)
				push!(TLR, (t1=df1.t[1],t2=df2.t[1],
						    x1=df1.x[1],y1=df1.y[1],z1=df1.z[1],
						    x2=df2.x[1],y2=df2.y[1],z2=df2.z[1]))

				#select reco
				hitdf  = JPetalo.reco_hits(event, ecut, pdf)

				# reconstruct (x,y) : barycenter
				b1, b2, hq2df, hq1df = fxy(hitdf)


				# find R1 and R2 (from True info)

				r1 = sqrt(df1.x[1]^2 + df1.y[1]^2)
				r2 = sqrt(df2.x[1]^2 + df2.y[1]^2)

				# New (x,y) positions estimated from R1 and R2
				x1,y1,z1  = JPetalo.radial_correction(b1, r1, rsipm)
				x2,y2,z2  = JPetalo.radial_correction(b2, r2, rsipm)

				push!(RLR, (t1=df1.t[1],t2=df2.t[1],
						    x1=x1,y1=y1,z1=z1,
						    x2=x2,y2=y2,z2=z2))

				push!(RLB, (t1=df1.t[1],t2=df2.t[1],
						    x1=b1.x,y1=b1.y,z1=b1.z,
						    x2=b2.x,y2=b2.y,z2=b2.z))

				cevt = event
			end
		end
	end
		return DataFrame(TLR[2:end]),
		       DataFrame(RLB[2:end]),
			   DataFrame(RLR[2:end])
end

"""
	df_to_mlemlor(ldf::DataFrame)

Take a LOR Data Frame and return a vector of MlemLor, suitable for hdf5
"""
function df_to_mlemlor(ldf::DataFrame)
	ml = [JPetalo.MlemLor(ldf.t1[i], ldf.t2[i],
		          ldf.x1[i],ldf.y1[i],ldf.z1[i],
				  ldf.x2[i],ldf.y2[i],ldf.z2[i]) for i in 1:nrow(ldf)]
	return ml
end

# """
# 	radial_correction(b::JPetalo.Hit, r::Float32, rsipm::Float32)
#
# Take the estimated radius of interaction (r), and the radius of the sipm
# and return the corrected positions
# """
# function radial_correction(b::JPetalo.Hit, r::Float32, rsipm::Float32)
# 	r2 = r / rsipm   # redial correction
# 	return r2 * b.x, r2 * b.y, b.z
# end

function makelors(args)
	ecut = 4.0

	println("+++makelor: ecut (in pes) = ", ecut)

	#dr = datadir("nema3-vac-1m")
	dr = datadir(args["dir"])
	files = glob("*.h5",dr)

	println("number of files in data dir = ", length(files))
	file_i = args["filei"]
	file_l = args["filel"]
	odir   = args["odir"]
	ecut   = args["ecut"]

	println("reading =", file_l - file_i + 1, " files")

	#tldf, rlhbdf, rlkbdf, rlkkdf = makelor_3p(files, 4.0, file_i, file_l)
	tldf, rlbdf, rlrdf = makelor_3p(files, JPetalo.lor_kmeans,
	                                ecut, file_i, file_l)

	mtl   = df_to_mlemlor(tldf)
	mrlb  = df_to_mlemlor(rlbdf)
	mrlr  = df_to_mlemlor(rlrdf)

	smtl  = string(odir,"/tl_",file_i,"_", file_l, ".h5")
	smrlb = string(odir,"/rlb_",file_i,"_", file_l, ".h5")
	smrlr = string(odir,"/rlr_",file_i,"_", file_l, ".h5")

	JPetalo.write_lors_hdf5(datadir(smtl), mtl)
	JPetalo.write_lors_hdf5(datadir(smrlb), mrlb)
	JPetalo.write_lors_hdf5(datadir(smrlr), mrlr)

end


function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--dir", "-d"
            help = "directory with nema3 simulations"
            arg_type = String
            default = "nema3-vac-1m"
		"--odir", "-o"
            help = "output directory"
            arg_type = String
            default = "nema3-vac-1m"

		"--filei", "-i"
	        help = "number of initial file in glob list"
	        default  = 1
			arg_type = Int
		"--filel", "-l"
		    help = "number of last file in glob list"
		    default  = 1
			arg_type = Int
		"--ecut", "-e"
		    help = "cut on SiPM charge (in pes)"
		    default  = 4.0
			arg_type = Float64

    end

    return parse_args(s)
end

function main()
	parsed_args = parse_commandline()
	println("Running makelors with arguments", parsed_args)
	makelors(parsed_args)
end

@time main()
