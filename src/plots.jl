using Plots
using LazySets

"""
	plot_truehits(th1, th2, emax=51.1)
"""
function plot_truehits(th1, th2, emax=51.1)

	pxyq = scatter(th1.x,th1.y,th1.e, legend=false)
	pxyq = scatter!(pxyq, th2.x,th2.y,th2.e, legend=false)


	pxy = scatter(th1.x,th1.y, marker_z = th1.e, markersize = th1.e/emax,
			         color = :jet, legend=false)
	pxy = scatter!(pxy, th2.x,th2.y, marker_z = th2.e, markersize = th2.e/emax,
			         color = :jet, legend=false)
	xlabel!("x")
	ylabel!("y")

	pxz = scatter(th1.x,th1.z, marker_z = th1.e/emax, markersize = th1.e/emax,
		         color = :jet,
		         legend=false)
	pxz = scatter!(pxz, th2.x,th2.z, marker_z = th2.e/emax, markersize = th2.e/emax,
		         color = :jet,
		         legend=false)
	xlabel!("x")
	ylabel!("z")

	pyz = scatter(th1.y,th1.z, marker_z = th1.e, markersize = th1.e/emax,
		         color = :jet,
		         legend=false)
	pyz = scatter!(pyz, th2.y,th2.z, marker_z = th2.e, markersize = th2.e/emax,
		         color = :jet,
		         legend=false)
	xlabel!("y")
	ylabel!("z")

	return pxyq,pxy,pxz,pyz
end


"""
	plot_lors(th1, th2, index, emax=51.1)
"""
function plot_lors(th1, th2, index, emax=51.1)
	h1 = select_truehit(th1, index)
	h2 = select_truehit(th2, index)


	pxy = scatter([h1.x,h2.x], [h1.y,h2.y], marker_z = [h1.e,h2.e],
		          markersize = [h1.e/emax, h2.e/emax],
			      color = :jet, legend=false)

	lsxy = LineSegment([h1.x,h1.y],[h2.x,h2.y])
	sxy = plot!(pxy,lsxy)
	xlabel!("x")
	ylabel!("y")

	pxz = scatter([h1.x,h2.x], [h1.z,h2.z], marker_z = [h1.e,h2.e],
		          markersize = [h1.e/emax, h2.e/emax],
			      color = :jet, legend=false)

	lsxz = LineSegment([h1.x,h1.z],[h2.x,h2.z])
	sxz = plot!(pxz,lsxz)
	xlabel!("x")
	ylabel!("z")

	pyz = scatter([h1.y,h2.y], [h1.z,h2.z], marker_z = [h1.e,h2.e],
		          markersize = [h1.e/emax, h2.e/emax],
			      color = :jet, legend=false)

	lsyz = LineSegment([h1.y,h1.z],[h2.y,h2.z])
	syz = plot!(pyz,lsyz)
	xlabel!("y")
	ylabel!("z")

	return sxy, sxz, syz
end


"""
	plot_lors_all(th1, th2, emax=51.1)
"""
function plot_lors_all2(th1, th2, emax=51.1)
	h1 = select_truehit(th1, 1)
	h2 = select_truehit(th2, 1)

	pxy  = scatter([h1.x,h2.x], [h1.y,h2.y], marker_z = [h1.e,h2.e],
		          markersize = [h1.e/emax, h2.e/emax],
			      color = :jet, legend=false)

	lsxy = LineSegment([h1.x,h1.y],[h2.x,h2.y])
	sxy  = plot!(pxy,lsxy)

	for indx in 2:size(th1.x)[1]
		h1 = select_truehit(th1, indx)
		h2 = select_truehit(th2, indx)

		pxy = scatter!(pxy, [h1.x,h2.x], [h1.y,h2.y], marker_z = [h1.e,h2.e],
			          markersize = [h1.e/emax, h2.e/emax],
				      color = :jet, legend=false)

	    lsxy = LineSegment([h1.x,h1.y],[h2.x,h2.y])
	    sxy = plot!(pxy,lsxy)
	end
	xlabel!("x")
	ylabel!("y")

	return sxy
end


function plot_lors_all(th1, th2, emax=100.0)
	function plot_xy()
		h1 = select_truehit(th1, 1)
		h2 = select_truehit(th2, 1)

		pxy  = scatter([h1.x,h2.x], [h1.y,h2.y], marker_z = [h1.e,h2.e],
					  markersize = [h1.e/emax, h2.e/emax],
					  color = :jet, legend=false)

		lsxy = LineSegment([h1.x,h1.y],[h2.x,h2.y])
		sxy  = plot!(pxy,lsxy)

		for indx in 2:size(th1.x)[1]
			h1 = select_truehit(th1, indx)
			h2 = select_truehit(th2, indx)

			pxy = scatter!(pxy, [h1.x,h2.x], [h1.y,h2.y], marker_z = [h1.e,h2.e],
						  markersize = [h1.e/emax, h2.e/emax],
						  color = :jet, legend=false)

			lsxy = LineSegment([h1.x,h1.y],[h2.x,h2.y])
			sxy = plot!(pxy,lsxy)
		end
		xlabel!("x")
		ylabel!("y")

		return sxy
	end
	function plot_xz()
		h1 = select_truehit(th1, 1)
		h2 = select_truehit(th2, 1)

		pxz  = scatter([h1.x,h2.x], [h1.z,h2.z], marker_z = [h1.e,h2.e],
					  markersize = [h1.e/emax, h2.e/emax],
					  color = :jet, legend=false)

		lsxz = LineSegment([h1.x,h1.z],[h2.x,h2.z])
		sxz  = plot!(pxz,lsxz)

		for indx in 2:size(th1.x)[1]
			h1 = select_truehit(th1, indx)
			h2 = select_truehit(th2, indx)

			pxz = scatter!(pxz, [h1.x,h2.x], [h1.z,h2.z], marker_z = [h1.e,h2.e],
						  markersize = [h1.e/emax, h2.e/emax],
						  color = :jet, legend=false)

			lsxz = LineSegment([h1.x,h1.z],[h2.x,h2.z])
			sxz = plot!(pxz,lsxz)
		end
		xlabel!("x")
		ylabel!("z")

		return sxz
	end

	function plot_yz()
		h1 = select_truehit(th1, 1)
		h2 = select_truehit(th2, 1)

		pyz  = scatter([h1.y,h2.y], [h1.z,h2.z], marker_z = [h1.e,h2.e],
					  markersize = [h1.e/emax, h2.e/emax],
					  color = :jet, legend=false)

		lsyz = LineSegment([h1.y,h1.z],[h2.y,h2.z])
		syz  = plot!(pyz,lsyz)

		for indx in 2:size(th1.x)[1]
			h1 = select_truehit(th1, indx)
			h2 = select_truehit(th2, indx)

			pyz = scatter!(pyz, [h1.y,h2.y], [h1.z,h2.z], marker_z = [h1.e,h2.e],
						  markersize = [h1.e/emax, h2.e/emax],
						  color = :jet, legend=false)

			lsyz = LineSegment([h1.y,h1.z],[h2.y,h2.z])
			syz = plot!(pyz,lsyz)
		end
		xlabel!("y")
		ylabel!("z")

		return syz
	end

	sxy = plot_xy()
	sxz = plot_xz()
	syz = plot_yz()
	return sxy, sxz, syz

end

"""
	plot_xyzq(hdf::DataFrame)

Take a hit dataframe and produce 4 scatter plots showing the clustering
of the SiPMs and the phi distribution showing how the clusters distribute
at opposite phi values
"""
function plot_xyzq(hdf, qnorm=50.0)
	pxyq = scatter(hdf.x,hdf.y,hdf.q, legend=false)
	pxy = scatter(hdf.x,hdf.y, marker_z = hdf.q, markersize = hdf.q/qnorm,
	              color = :jet,
		          legend=false)
	xlabel!("x")
	ylabel!("y")
	pxz = scatter(hdf.x,hdf.z, marker_z = hdf.q, markersize = hdf.q/qnorm,
	             color = :jet,
		         legend=false)
	xlabel!("x")
	ylabel!("z")
	pyz = scatter(hdf.y,hdf.z, marker_z = hdf.q, markersize = hdf.q/qnorm,
	              color = :jet,
		          legend=false)
	xlabel!("y")
	ylabel!("z")

	return pxyq, pxy, pxz, pyz
end

"""
	plot_barycenter(bp,bn, pxy, pxz, pyz)
"""
function plot_barycenter(bp, bn, pxy, pxz, pyz, qnorm=1.0)
	lsxy = LineSegment([bp.x,bp.y],[bn.x,bn.y])
	sxy = scatter!(pxy,[bp.x,bn.x],[bp.y,bn.y],
		marker_z = [bp.q,bn.q], markersize = qnorm,  color = :jet, legend=false)
	sxy = plot!(sxy,lsxy)
	syz = scatter!(pyz,[bp.y,bn.y],[bp.z,bn.z],
		marker_z = [bp.q,bn.q], markersize = qnorm, color = :jet,legend=false)
	lsyz = LineSegment([bp.y,bp.z],[bn.y,bn.z])
	syz  = plot!(syz,lsyz)
	sxz = scatter!(pxz,[bp.x,bn.x],[bp.z,bn.z],
		marker_z = [bp.q,bn.q], markersize = qnorm, color = :jet,legend=false)
	lsxz = LineSegment([bp.x,bp.z],[bn.x,bn.z])
	sxz  = plot!(sxz,lsxz)
	return sxy, syz, sxz
end


function plot_lors_barycenter(BP, BN, qmax=51.1)

	function plot_xy(BP,BN)
		h1 = BP[1]
		h2 = BN[1]

		pxy  = scatter([h1.x,h2.x], [h1.y,h2.y], marker_z = [h1.q,h2.q],
		               markersize = [h1.q/qmax, h2.q/qmax],
			           color = :jet, legend=false)

		lsxy = LineSegment([h1.x,h1.y],[h2.x,h2.y])
		sxy  = plot!(pxy,lsxy)

		for i in 2:size(BP)[1]
			h1 = BP[i]
			h2 = BN[i]

			pxy = scatter!(pxy, [h1.x,h2.x], [h1.y,h2.y], marker_z = [h1.q,h2.q],
		     	          markersize = [h1.q/qmax, h2.q/qmax],
			     	      color = :jet, legend=false)

	    	lsxy = LineSegment([h1.x,h1.y],[h2.x,h2.y])

	    	sxy = plot!(pxy,lsxy)
		end
		return sxy
	end

	function plot_xz(BP,BN)
		h1 = BP[1]
		h2 = BN[1]

		pxz  = scatter([h1.x,h2.x], [h1.z,h2.z], marker_z = [h1.q,h2.q],
		               markersize = [h1.q/qmax, h2.q/qmax],
			           color = :jet, legend=false)

		lsxz = LineSegment([h1.x,h1.z],[h2.x,h2.z])
		sxz  = plot!(pxz,lsxz)

		for i in 2:size(BP)[1]
			h1 = BP[i]
			h2 = BN[i]

			pxz = scatter!(pxz, [h1.x,h2.x], [h1.z,h2.z], marker_z = [h1.q,h2.q],
		     	          markersize = [h1.q/qmax, h2.q/qmax],
			     	      color = :jet, legend=false)

	    	lsxz = LineSegment([h1.x,h1.z],[h2.x,h2.z])

	    	sxz = plot!(pxz,lsxz)
		end
		return sxz
	end

	function plot_yz(BP,BN)
		h1 = BP[1]
		h2 = BN[1]

		pyz  = scatter([h1.y,h2.y], [h1.z,h2.z], marker_z = [h1.q,h2.q],
		               markersize = [h1.q/qmax, h2.q/qmax],
			           color = :jet, legend=false)

		lsyz = LineSegment([h1.y,h1.z],[h2.y,h2.z])
		syz  = plot!(pyz,lsyz)

		for i in 2:size(BP)[1]
			h1 = BP[i]
			h2 = BN[i]

			pyz = scatter!(pyz, [h1.y,h2.y], [h1.z,h2.z], marker_z = [h1.q,h2.q],
		     	          markersize = [h1.q/qmax, h2.q/qmax],
			     	      color = :jet, legend=false)

	    	lsyz = LineSegment([h1.y,h1.z],[h2.y,h2.z])

	    	syz = plot!(pyz,lsyz)
		end
		return syz
	end


	sxy = plot_xy(BP,BN)
	xlabel!("x")
	ylabel!("y")

	sxz = plot_xz(BP,BN)
	xlabel!("x")
	ylabel!("z")

	syz = plot_yz(BP,BN)
	xlabel!("y")
	ylabel!("z")

	return sxy, sxz, syz
end
