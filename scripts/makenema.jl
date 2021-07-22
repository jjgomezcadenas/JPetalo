push!(LOAD_PATH,"../src/")
using DrWatson
@quickactivate "JPetalo"
include("../src/jpetalo.jl")

using DataFrames
using CSV
using Glob
using ArgParse
using Logging

logger = SimpleLogger(stdout, Logging.Warn)
old_logger = global_logger(logger)


function makenema(args)

	qc = Float32(args["qcut"])
	pde = Float32(args["pde"])
	qmin = Float32(args["qmin"])
	qmax = Float32(args["qmax"])
	cq1  = Float32(args["cq1"])
	cq2  = Float32(args["cq2"])
	maxpes = args["maxpes"]
	ntof = args["ntof"]
	sigmatof = Float32(args["sigmatof"])
	lorf = args["loralgo"]

	dr = datadir(args["dir"])
	outd = datadir(args["odir"])
	outf = args["ofile"]
	file_i = args["filei"]
	file_l = args["filel"]

	lor_algo = JPetalo.lor_kmeans

	if lorf == "lor_qmax"
		lor_algo = JPetalo.lor_maxq
	end

	output = string(outd,"/", outf)
	files = glob("*.h5",dr)
	phot =args["phot"]

	@info "makenema configuration"
	@info "pde  = $pde sigma tof (ps) = $sigmatof maxpes=$maxpes"
	@info "ecut (pes)  = $qc qmin (pes) = $qmin qmx (pes)=$qmax"
	@info " r = f(q) parameters: cq1  = $cq1 cq2 = $cq2"
	@info " lor function cq1  = $lorf"
	@info " photoelectric only  = $phot"


	@info("number of files in data dir = $length(files)")
	@info("reading = $(file_l - file_i + 1) files")
	@info("output file  = $output")

	n3df = JPetalo.nemareco(files, file_i, file_l,
	                        pde, maxpes, sigmatof,
	                        qc, qmin, qmax, cq1, cq2, ntof, phot, lor_algo)

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
			arg_type = Int
			default = 10
		"--ntof", "-c"
			help = "number of sipms for average"
			arg_type = Int
			default = 5
		"--prteach", "-e"
			help = "print each events"
			arg_type = Int
			default = 1000
		"--sigmatof", "-t"
			help = "smearing on TOF due to electronics and sensor"
			arg_type = Float64
			default = 0.085  # in ns
		"--loralgo", "-g"
			help = "algorithm to use for LOR reconstruction "
			arg_type = String
			default = "lor_kmeans"
		"--qmin", "-m"
			help = "qmin on SiPMs"
			arg_type = Float64
			default = 1800.
		"--qmax", "-M"
            help = "qmax on SiPMs"
			arg_type = Float64
            default = 3000.
		"--cq1", "-I"
            help = "intercept in function r =cq1 + cq2*q"
			arg_type = Float64
            default = 297.9
		"--cq2", "-S"
            help = "slope in function r =cq1 + cq2*q"
			arg_type = Float64
            default = 0.031
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
	makenema(parsed_args)
end

@time main()
