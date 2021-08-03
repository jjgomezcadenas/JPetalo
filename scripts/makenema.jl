push!(LOAD_PATH,"../src/")
using DrWatson
#@quickactivate(@__DIR__)
#@quickactivate "JPetalo"
quickactivate("../.", "JPetalo")
include("../src/JPetalo.jl")

using DataFrames
using CSV
using Glob
using ArgParse
using Logging
using Printf

logger = SimpleLogger(stdout, Logging.Warn)
old_logger = global_logger(logger)

# radius-from-q function. The parameters of the straight line have been
# obtained from a fit to photoelectric data.

# function grfq(p::Vector{Float32})
# 	function rfq(q)
# 		return p[1] + p[2] * q
# 	end
# 	return rfq
# end
#
# # radius-from-zstd function. The parameters of the parabol have been
# # obtained from a fit to photoelectric data.
#
# function grfz(p::Vector{Float32})
# 	function rfz(z)
# 		return p[1] + p[2] * z + p[3] * z^2
# 	end
# 	return rfz
# end

function makenema(args)

	lorf    = args["loralgo"]
	detconf = args["detconf"]
	phot    = args["phot"]
	dr      = datadir(args["dir"])
	outd    = datadir(args["odir"])
	outf    = args["ofile"]
	file_i  = args["filei"]
	file_l  = args["filel"]

	lor_algo = JPetalo.lor_kmeans
	if lorf == "lor_qmax"
		lor_algo = JPetalo.lor_maxq
	end

	output = string(outd,"/", outf)
	files = glob("*.h5",dr)

	if detconf == "pde_1_sigmatof_1ps"
		pde  = 1.0f0
		sigma_tof = 0.001f0
		ecut = 10.0f0
		qmin = 100.0f0
		qmax = 15000.0f0
		max_pes = 10
		ntof =5
		#rfq = grfq([315.1f0, 0.008f0])
		#rfz = grfq([392.5f0, 0.36f0, -0.065f0])
		#dconf = JPetalo.DetConf(pde, sigma_tof, ecut, qmin, qmax, max_pes, ntof, rfq, rfz)
		dconf = JPetalo.DetConf(pde, sigma_tof, ecut, qmin, qmax, max_pes, ntof, rfq, rfz)
	elseif detconf == "pde_0.3_sigmatof_1ps"
		pde  = 0.3f0
		sigma_tof = 0.001f0
		ecut = 3.0f0
		qmin = 100.0f0
		qmax = 5000.0f0
		max_pes = 10
		ntof =5
		dconf = JPetalo.DetConf(pde, sigma_tof, ecut, qmin, qmax, max_pes, ntof)
	elseif detconf == "pde_0.3_sigmatof_85ps"
		pde  = 0.3f0
		sigma_tof = 0.085f0
		ecut = 3.0f0
		qmin = 100.0f0
		qmax = 5000.0f0
		max_pes = 10
		ntof =5
		dconf = JPetalo.DetConf(pde, sigma_tof, ecut, qmin, qmax, max_pes, ntof)
	elseif detconf == "pde_0.3_sigmatof_85ps_qcuts"
		pde  = 0.3f0
		sigma_tof = 0.085f0
		ecut = 3.0f0
		qmin = 1500.0f0
		qmax = 3000.0f0
		max_pes = 10
		ntof =5
		dconf = JPetalo.DetConf(pde, sigma_tof, ecut, qmin, qmax, max_pes, ntof)
	else                       # by default, 40 ps jitter  & 30 ps electronics = 50 ps
		pde  = 0.3f0
		sigma_tof = 0.050f0
		ecut = 5.0f0           # increase ecut 
		qmin = 1500.0f0
		qmax = 3000.0f0
		max_pes = 10
		ntof =7                 # increase sipms for average
		dconf = JPetalo.DetConf(pde, sigma_tof, ecut, qmin, qmax, max_pes, ntof)

	end

	println("makenema configuration")
	println("data in", datadir())
	println("detector configuration", dconf)
	println(" photoelectric only  = $phot")
	println(" lor_algo  = $lorf")

	println("number of files in data dir = $(length(files))")
	println("reading = $(file_l - file_i + 1) files")
	println("output file  = $output")

	n3df = JPetalo.nemareco(files, dconf, file_i, file_l, phot, lor_algo)

	CSV.write(output, n3df)
end


function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--dir", "-d"
            help = "directory with nema simulations"
            arg_type = String
            default = "nema3-vac-1m"
		"--odir", "-o"
            help = "output directory"
            arg_type = String
            default = "evtdf"
		"--ofile", "-x"
            help = "output file"
            arg_type = String
            default = "evtdf.csv"
		"--filei", "-i"
	        help = "number of initial file in glob list"
	        default  = 1
			arg_type = Int
		"--filel", "-l"
		    help = "number of last file in glob list"
		    default  = 1
			arg_type = Int
		"--loralgo", "-g"
			help = "algorithm to use for LOR reconstruction "
			arg_type = String
			default = "lor_kmeans"
		"--detconf", "-c"
			help = "Detector configuration"
			arg_type = String
			default = "all"
		"--phot", "-p"
			help = "Select photoelectric if 1"
			action = :store_true
    end

    return parse_args(s)
end

function main()
	parsed_args = parse_commandline()
	println("Running makenema with arguments", parsed_args)
	makenema(parsed_args)
end

@time main()
