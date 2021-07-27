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



function makelors(args)

	dr = datadir(args["dir"])
	files = glob("*.csv",dr)

	@info "directory contains $length(files))"

	file_i = args["filei"]
	file_l = args["filel"]
	odir   = args["odir"]

	for file in files[file_i:file_l]               # loop on files
		@info "reading file = $file"

		n3df  = DataFrame(CSV.File(file))
		n3dfu = setunits(n3df)



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
            help = "directory with nema dfs"
            arg_type = String
            default = "nema3all"
		"--odir", "-o"
            help = "output directory"
            arg_type = String
            default = "lorsall"

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
	println("Running makelors with arguments", parsed_args)
	makelors(parsed_args)
end

@time main()
