push!(LOAD_PATH,"../src/")
using DrWatson
@quickactivate "JPetalo"
include("../src/jpetalo.jl")

using DataFrames
using CSV
using Glob
using ArgParse


function makenema3(args)

	qc = Float32(args["qcut"])
	pde = Float32(args["pde"])
	qmin = Float32(args["qmin"])
	qmax = Float32(args["qmax"])
	dr = datadir(args["dir"])
	outd = datadir(args["odir"])
	outf = args["ofile"]
	file_i = args["filei"]
	file_l = args["filel"]

	output = string(outd,"/", outf)
	files = glob("*.h5",dr)

	phot =args["phot"]
	warn =args["warn"]


	println("number of files in data dir = ", length(files))
	println("reading =", file_l - file_i + 1, " files")
	println("output file  = ", output)
	println("ecut (in pes) = ", qc)
	println("pde  = ", pde)
	println("photoelectric only  = ", phot)
	println("warn  = ", warn)


	n3df = JPetalo.nema3a(files, file_i, file_l, qc, pde, qmin, qmax, phot, warn)

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
		"--warn"
	 		help = "activate warnings"
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
