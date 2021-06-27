using DataFrames
using Clustering
using StatsModels

"""
	reco_lor(pdf::petaloDF, evtlist::Vector{Int64}, ecut::Number)

	Computes reconstructed lors

"""
function reco_lor(pdf::PetaloDF, evtlist::Vector{Integer}, ecut::Float64)

	function ksipmsel(hdf::DataFrame, ka::Vector{Int64})
		return hdf[(ka.==2), :], hdf[(ka.==1), :]
	end

	function get_hits_as_matrix(hitdf)
		f = @formula(0 ~ x + y + z)
		f2 = apply_schema(f, schema(f, hitdf))
		resp, pred = modelcols(f2, hitdf)
		return transpose(pred)
	end

	BP  = []
	BN  = []
	rBP = []
	rBN = []
	QM  = []
	for evt0 in evtlist
		evt    = select_event(pdf.total_charge, evt0)    #select event
		evtQ1  = evt[evt.charge.>ecut,:]                 # charge cut
		hitQdf = sipm_xyzq(evtQ1, pdf.sensor_xyz)        # reco hits
		qm     = maximum(hitQdf.q)

		# kmeans algo
		Mhits = get_hits_as_matrix(hitQdf)
		kr = kmeans(Mhits, 2)
		ka = assignments(kr) # get the assignments of points to clusters
		kc = counts(kr) # get the cluster sizes
		rk = kr.centers # get the cluster centers

		hq2df, hq1df = ksipmsel(hitQdf, ka)   # select using kmeans list
		b1 = baricenter(hq1df)     # baricenters
		b2 = baricenter(hq2df)

		rb1 = Hit(rk[1,1],rk[2,1], rk[3,1], b1.q)  #cluster centers
		rb2 = Hit(rk[1,2],rk[2,2], rk[3,2], b2.q)

		push!(BP, b1)
		push!(BN, b2)
		push!(rBP, rb1)
		push!(rBN, rb2)
		push!(QM, qm)

	end
	return QM, BP, BN, rBP, rBN
end
