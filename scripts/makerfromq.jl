push!(LOAD_PATH,"../src/")
using DrWatson
@quickactivate "JPetalo"
include("../src/jpetalo.jl")

using DataFrames
using CSV
using Glob
using ArgParse

function grfq(q)
    return 326.9 + 0.0226 * q  # replace these parameters by cqr[1], cqr[2]
end

rfq(q::Float32) = grfq(q)
rfq(q::Float64) = grfq(q)

function makerfq(args)

	dr = datadir(args["dir"])
	fl = args["file"]
	qmin = args["qmin"]
	qmax = args["qmax"]
	nbin = args["nbin"]
	println("qmin = ", qmin, " qmax = ", qmax)

	input = string(dr,"/",fl)
	n3df = DataFrame(CSV.File(input))

	# select events in the interval (qmin, qmax)
	n3dfq = JPetalo.select_by_column_value_interval(n3df, "q1", qmin, qmax)
	pqrdf = JPetalo.p1df(n3dfq.q1, n3dfq.r, nbin)
	lfqr, pqr, cqr = JPetalo.lfit(pqrdf)
	println("parameters from linear fit (r  = c1 + c2 * q)")
	println("c1 =", cqr[1], " c2 =", cqr[2])

end


function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--dir", "-d"
            help = "directory with nema3df data"
            arg_type = String
            default = "nema3df"
		"--file", "-f"
            help = "file with nema3df data"
            arg_type = String
            default = "nema3df_f600_q4_1300_3000_phot.csv"
		"--qmin", "-m"
			help = "qmin on SiPMs"
			arg_type = Float64
			default = 1400.
		"--qmax", "-M"
            help = "qmax on SiPMs"
			arg_type = Float64
            default = 3000.
		"--nbin", "-n"
            help = "number of bins for profile fit"
			arg_type = Int
            default = 100
    end

    return parse_args(s)
end

function main()
	parsed_args = parse_commandline()
	println("Running makerfq with arguments", parsed_args)
	makerfq(parsed_args)
end

@time main()
