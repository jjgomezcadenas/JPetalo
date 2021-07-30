push!(LOAD_PATH,"../src/")
using DrWatson
@quickactivate "JPetalo"
include("../src/JPetalo.jl")

import .JPetalo: Dtsel, dtfirst, dtminimum, dtaverage, postrue, posreco
import .JPetalo: Possel, setunits, dftolor, write_lors_hdf5

using DataFrames
using CSV
using Glob
using ArgParse
using Logging

logger = SimpleLogger(stdout, Logging.Warn)
old_logger = global_logger(logger)

dtselmap =Dict("first"=>dtfirst, "minimum"=>dtminimum, "average"=>dtaverage)
posselmap =Dict("true"=>postrue, "reco"=>posreco)

function makelors(args)

	dr = datadir(args["dir"])
	files = glob("*.csv",dr)

	@info "directory contains $length(files))"

	file_i = args["filei"]
	file_l = args["filel"]
	odir   = args["odir"]
	outf   = args["ofile"]
	selt   = dtselmap[args["selt"]]
	selp   = posselmap[args["selpos"]]

	output = string(odir,"/", outf, "-", args["selt"], "-", args["selpos"], ".h5")
	println("writing output to $output")

	LORS = []
	for file in files[file_i:file_l]               # loop on files
		@info "reading file = $file"
		evtdfu = setunits(DataFrame(CSV.File(file)))
		lors   = dftolor(evtdfu, selt, selp)
		push!(LORS, lors)
	end

	write_lors_hdf5(datadir(output), reduce(vcat, LORS))

end

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--selt", "-t"
            help = "select time to be used in lor (dtfirst, dtmimimum, dtaverage)"
            arg_type = String
            default = "first"
        "--selpos", "-p"
            help = "select position to be used in lor (postrue, posreco)"
            arg_type = String
            default = "true"
        "--dir", "-d"
            help = "directory with nema dfs"
            arg_type = String
            default = "evtdfall"
		"--odir", "-o"
            help = "output directory"
            arg_type = String
            default = "lors"
		"--filei", "-i"
	        help = "number of initial file in glob list"
	        default  = 1
			arg_type = Int
		"--filel", "-l"
		    help = "number of last file in glob list"
		    default  = 1
			arg_type = Int
		"--ofile", "-x"
            help = "output file prefix"
            arg_type = String
            default = "all"
    end

    return parse_args(s)
end

function main()
	parsed_args = parse_commandline()
	println("Running makelors with arguments", parsed_args)
	makelors(parsed_args)
end

@time main()
