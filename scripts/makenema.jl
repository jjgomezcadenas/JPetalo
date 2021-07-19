push!(LOAD_PATH,"../src/")
using DrWatson
@quickactivate "JPetalo"
include("../src/jpetalo.jl")

using DataFrames
using CSV
using Glob
using ArgParse
using Logging

logger = SimpleLogger(stdout, Logging.Debug)
old_logger = global_logger(logger)


function makenema(args)

	qc = Float32(args["qcut"])
	pde = Float32(args["pde"])
	qmin = Float32(args["qmin"])
	qmax = Float32(args["qmax"])
	maxpes = args["maxpes"]
	prteach = args["prteach"]
	sigmatof = Float32(args["sigmatof"])
	lorf = args["loralgo"]
	dr = datadir(args["dir"])
	outd = datadir(args["odir"])
	outf = args["ofile"]
	file_i = args["filei"]
	file_l = args["filel"]

	lor_algo = JPetalo.lor_maxq

	if lorf == "lor_kmeans"
		lor_algo = JPetalo.lor_kmeans
	end

	output = string(outd,"/", outf)
	files = glob("*.h5",dr)
	phot =args["phot"]

	println(" makenema configuration")
	println("pde  = ", pde)
	println("sigma tof  = ", sigmatof)
	println("ecut (in pes) = ", qc)
	println("qmin (total charge in pes) = ", qmin)
	println("qmax (total charge in pes) = ", qmax)
	println("lor function  = ", lorf)
	println("photoelectric only  = ", phot)

	println("number of files in data dir = ", length(files))
	println("reading =", file_l - file_i + 1, " files")
	println("output file  = ", output)

	n3df = JPetalo.nema3a(files, file_i, file_l,
	                      pde, maxpes, sigmatof,
	                      qc, qmin, qmax, prteach, phot, lor_algo)


	CSV.write(output, n3df)
end


function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--dir", "-d"
            help = "directory with nema3 simulations"
            arg_type = String
            default = "nema3-vac-1m"
		"--qcut", "-q"
            help = "cut on q of SiPMs"
			arg_type = Float64
            default = 2.0
		"--pde", "-p"
            help = "PDE of SiPMs"
			arg_type = Float64
            default = 0.3
		"--maxpes", "-s"
			help = "max number of pes stores in a SiPM in the waveform"
			arg_type = Integer
			default = 10
		"--prteach", "-e"
			help = "print each events"
			arg_type = Integer
			default = 1000
		"--sigmatof", "-t"
			help = "smearing on TOF due to electronics and sensor"
			arg_type = Float64
			default = 0.085  # in ns
		"--loralgo", "-l"
			help = "algorithm to use for LOR reconstruction "
			arg_type = String
			default = "lor_maxq"
		"--qmin", "-m"
			help = "qmin on SiPMs"
			arg_type = Float64
			default = 1400.
		"--qmax", "-M"
            help = "qmax on SiPMs"
			arg_type = Float64
            default = 3000.
		"--phot"
	 		help = "Select photoelectric if 1"
	 		action = :store_true
		"--odir", "-o"
            help = "output directory"
            arg_type = String
            default = "nema3df"
		"--ofile", "-x"
            help = "output file"
            arg_type = String
            default = "nema3df.csv"
		"--filei", "-i"
	        help = "number of initial file in glob list"
	        default  = 1
			arg_type = Int
		"--filel", "-l"
		    help = "number of last file in glob list"
		    default  = 1
			arg_type = Int
    end

    return parse_args(s)
end

function main()
	parsed_args = parse_commandline()
	println("Running makenema with arguments", parsed_args)
	makenema3(parsed_args)
end

@time main()
