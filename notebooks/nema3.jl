### A Pluto.jl notebook ###
# v0.14.5

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

# ╔═╡ d16a3879-acb4-4cce-9b87-faadf4abfece
import Pkg; Pkg.add.(["HTTP"])

# ╔═╡ 4ee98da9-f3ce-4782-9538-f878f27ed9f7
Pkg.add.(["HDF5"])

# ╔═╡ 2e5dd476-4e23-427e-8e72-a6af044eb397
Pkg.add.(["Clustering"])

# ╔═╡ 26a0d201-f36c-4bcf-9219-a2f22c500598
Pkg.add.(["StatsModels"])

# ╔═╡ 3207f446-8643-11eb-37ba-c9aec47fcb8f
begin
	using Markdown
	using InteractiveUtils
	using PlutoUI
	using Test
	using Plots
	using StatsPlots
	using QuadGK
	using CSV
	using DataFrames
	using Statistics
	using StatsPlots
	using LinearAlgebra
	using HTTP
	using HDF5
	#import Pandas
end

# ╔═╡ 5115917a-8644-11eb-19fc-0528741ca75d
begin
	using Unitful
	using UnitfulEquivalences
	using PhysicalConstants.CODATA2018
	using LazySets
	using Formatting
	using Printf
	using StrLiterals
	using StrFormat
	using Format
end

# ╔═╡ fc8b79a2-8728-11eb-2da7-e3ffa3ceef08
using DrWatson

# ╔═╡ ac019a72-dca3-4e7e-9ffe-3cfbede4106e
using Clustering

# ╔═╡ aa9f1b97-21cb-4592-8f47-8549f1305da2
using StatsModels

# ╔═╡ 79cfd2fc-9046-11eb-2b13-1b877d57645d
md"# NEMA3

- NEMA3 studies
"

# ╔═╡ 0068d015-f292-4bf7-82f9-c6c0115f96e2
#using Plotly

# ╔═╡ 68e738e7-88bd-41c2-89e4-594f07d64ddc
function ingredients(path::String)
	# this is from the Julia source code (evalfile in base/loading.jl)
	# but with the modification that it returns the module instead of the last object
	name = Symbol(basename(path))
	m = Module(name)
	Core.eval(m,
        Expr(:toplevel,
             :(eval(x) = $(Expr(:core, :eval))($name, x)),
             :(include(x) = $(Expr(:top, :include))($name, x)),
             :(include(mapexpr::Function, x) = $(Expr(:top, :include))(mapexpr, $name, x)),
             :(include($path))))
	m
end

# ╔═╡ 0f2f4c78-8729-11eb-2bab-27812ce8c47e
@quickactivate "JPetalo"

# ╔═╡ 5ee27d52-86fd-11eb-365e-9f2b1a095575
;pwd()

# ╔═╡ 621ec96c-86fd-11eb-1c41-379cc17180dc
datadir()

# ╔═╡ 9b853f27-4288-42a6-8f12-ca004e1773b7
srcdir()

# ╔═╡ 26225a47-c7aa-4319-89f2-1abbc5c1116a
function void()
	println("")
end

# ╔═╡ cf89b973-7b0f-483d-8ce9-ba426f1df2a6
jp = ingredients(srcdir("jpetalo.jl"))

# ╔═╡ 0b1cedc7-8ada-45a5-adef-fbae794dee3e
markercolors = [:green :orange :black :purple :red  :yellow :brown :white]

# ╔═╡ 5fb08873-1ca6-44f2-b68a-438fea6007ed
gr(size=(500,500), xtickfontsize=10, ytickfontsize=10, xguidefontsize=10, yguidefontsize=10, legendfontsize=10, dpi=100, grid=(:y, :gray, :solid, 1, 1.0));


# ╔═╡ 383048a4-12e4-492f-bcbd-a83c6e7dc7aa
Plots.GRBackend()

# ╔═╡ 80542fd1-843e-4d78-9bd3-169c1d6a9672
md"# Notebook"

# ╔═╡ edbcd276-6e0d-41f3-91d8-131b0fc7486b
md"### Read NEMA3 DST"

# ╔═╡ fdff3aeb-834c-4726-8e6c-a35e6b9d46d1
path = datadir("nema3/nema3-vac-25k-1.h5")

# ╔═╡ bc16f5c1-863e-47de-87aa-838578157188
fid = h5open(path, "r") 

# ╔═╡ 5b8468a8-4465-4708-afd0-9c7495b3a8a3
pdf = jp.JPetalo.read_abc(path)

# ╔═╡ 8db332bb-9a40-4bdf-b018-64f6a35f6bdf
md"- ### Position of the SiPMs

- The scatter plot shows that the positions of the SiPMs form a ring in XY
- And a wall in XZ
" 

# ╔═╡ b8aacd81-9d8f-4105-a214-220be62b3e43
psipm_xy = scatter(pdf.sensor_xyz.x, pdf.sensor_xyz.y, title = "XY positions SiPMs", leg=false, fmt = :png);

# ╔═╡ 42e1caa1-4aba-44e7-ab32-4538662bcf29
psipm_xz = scatter(pdf.sensor_xyz.x, pdf.sensor_xyz.z, title = "XZ positions SiPMs", leg=false, fmt = :png);

# ╔═╡ fc41144a-ba14-4543-9523-c5450744e126
begin
gr(size=(500,500))
plot(psipm_xy,psipm_xz, layout = (1, 2), aspect_ratio=:equal, xlim=(-500.0, 500.0), 
		ylim=(-500.0, 500.0), size = (1400, 1000), legend=false, fmt = :png)
end

# ╔═╡ a2f23d32-150b-4bb7-96d0-79f687747e40
md" - Direction cosines of primary gammas"

# ╔═╡ 6fd444d3-03f3-4b56-a02d-cddba85ca96f
begin
	pvert_xy = Plots.histogram2d(pdf.primaries.vx, pdf.primaries.vy, 
		               title = "VXY ", 	
		               leg=false, fmt = :png);
	pvert_xz = Plots.histogram2d(pdf.primaries.vx, pdf.primaries.vz, title = "XZ ", 
		               leg=false, fmt = :png);
	pvert_yz = Plots.histogram2d(pdf.primaries.vy, pdf.primaries.vz, title = "YZ ", 
		               leg=false, fmt = :png);
	true
end

# ╔═╡ 06208989-8940-4cd5-8f59-c5e7c167b97b
plot(pvert_xy,pvert_xz, pvert_yz, layout = (1, 3), aspect_ratio=:equal,
	xlim=(-0.5, 0.5), ylim=(-1.0, 1.0), size = (1400, 1000), 
	legend=false, fmt = :png);

# ╔═╡ 5c4b2813-991e-4908-833c-195309dd95f1
md" - Position of sources"

# ╔═╡ bfc9e52f-7889-48fc-a6a8-55b8844420e0
begin
	p_xy = Plots.scatter(pdf.primaries.x, pdf.primaries.y, title = "XY primary ", 	
		               leg=false, fmt = :png);
	p_xz = Plots.scatter(pdf.primaries.x, pdf.primaries.z, title = "XZ primary ", 
		               leg=false, fmt = :png);
	p_yz = Plots.scatter(pdf.primaries.y, pdf.primaries.z, title = "YZ primary ", 
		               leg=false, fmt = :png);
	true
end

# ╔═╡ abfb7e89-99dd-4e2a-8b5f-ec8fefb8fdc3
plot(p_xy,p_xz, p_yz, layout = (1, 3), aspect_ratio=:equal,size = (1400, 1000), legend=false, fmt = :png)

# ╔═╡ 766c0b57-9a16-44d1-ab5f-e7c1d39a5176
md"- Each data file is indexed by event\_id and sensor\_id. The field sensor_id is used to extract the positions of the sipms in the DB and also provides the charge for a given time bin. This analysis is not considering time (in fact the data is binned with very large data bins, which in practice are irrelevant). Time analysis is illustrated in another notebook"

# ╔═╡ 48a21c9c-cfc3-4d23-b1d9-cb9e47753bb4
md"# Analysis"

# ╔═╡ 361424ff-0d5d-4150-8928-928a184db66f
md"### True info"

# ╔═╡ 2231ac80-f3cf-446b-b6dd-ba02efbd29b7
md"- Select events in LXe"

# ╔═╡ 9e792d68-79a7-4418-b96f-13c32be0e20c
vlxe = jp.JPetalo.select_by_column_value(pdf.vertices, "volume_id", 0);

# ╔═╡ 868a8593-f112-43ed-af79-88d972fd2791
md"- Select photoelectric"

# ╔═╡ 3e514c2f-03b9-42da-8317-0c80309066a9
vlxephe = jp.JPetalo.select_by_column_value(vlxe, "process_id", 1);

# ╔═╡ e025a1d1-d8ef-4ee7-8c7a-a58eb4d6eee2
md"- Select parent_id = 0 (primaries)"

# ╔═╡ cbd4bcba-1429-4bcf-8fec-20316b0f18e2
vlxephepr = jp.JPetalo.select_by_column_value(vlxephe, "parent_id", 0);

# ╔═╡ b8d08354-9346-4e70-ab23-9a7927524099
md"- Select by event (example)"

# ╔═╡ 0a7cfeff-e1e4-4662-bc59-5088e95749b6
df2 = jp.JPetalo.select_by_column_value(vlxephepr, "event_id", 10041);

# ╔═╡ 6c355936-7d4a-456c-bddb-0de22d3144fe
GL = jp.JPetalo.true_lors(vlxephepr);

# ╔═╡ bb427097-45b5-42a3-9c4b-ab48913187bd
th1, th2 = jp.JPetalo.get_truehits(GL);

# ╔═╡ ce932b08-2e26-4527-94f1-da3c49a5433b
th1.event_id

# ╔═╡ 99500b8e-3559-495e-83e9-42dde08173b6
th2.event_id

# ╔═╡ 9963d767-a54b-40db-ba6a-a3d64f6dcbd7
pxyqt,pxyt,pxzt,pyzt = jp.JPetalo.plot_truehits(th1, th2, 51.1);

# ╔═╡ df6d7369-88ae-4536-b66c-369d56291bce
plot(pxyqt)

# ╔═╡ 21d2c293-ec33-48a7-83f4-8cd6df8c0a70
plot(pxyt,pxzt,pyzt, layout = (1, 3), aspect_ratio=:equal,size = (1400, 1000), legend=false, fmt = :png)

# ╔═╡ 4d24d7fa-28f7-4d0a-9889-f0790ae6fe53
tlxy, tlxz, tlyz = jp.JPetalo.plot_lors(th1, th2, 1);

# ╔═╡ 302b25f2-cd53-4918-983d-dc37399507a3
plot(tlxy, tlxz, tlyz, layout = (1, 3), aspect_ratio=:equal,size = (1400, 1000), legend=false, fmt = :png)

# ╔═╡ 178d5a6c-2fb8-49cc-884b-414f86ceaad6
ptsxy, ptsxz, ptsyz = jp.JPetalo.plot_lors_all(th1, th2,101.)

# ╔═╡ 2694fc3f-2883-4af7-ba40-6a5b41ce1147
plot(ptsxy, layout= (1, 1), aspect_ratio=:equal,size = (1400, 1000),legend = false,  fmt = :png)

# ╔═╡ f7b06fae-650c-4ea2-9a12-4a28b775f845
plot(ptsxz, layout= (1, 1), size = (1600, 800),
	ylims=(-100.0, 100.0), legend = false,  fmt = :png)

# ╔═╡ f1d2767f-c729-4f22-8e83-5cd3161f7d93
plot(ptsyz, layout= (1, 1), size = (1600, 800),
	ylims=(-100.0, 100.0), legend = false,  fmt = :png)

# ╔═╡ 8b050930-1c99-4613-aa95-90f40743acab


# ╔═╡ f4f30742-121d-49ea-ad14-df9650bf6c6b
md"### Reconstructed data"

# ╔═╡ 801ded29-0d7f-4407-94e8-7fffb95613b8
th1

# ╔═╡ 8a6c603d-a33b-440b-b5bd-9b7bef37b2b5
md" - First step: select an event"

# ╔═╡ f3b71ecf-fa44-480d-bb7d-612edc6ff15e
@bind evt0 NumberField(1000:50000; default=10038)

# ╔═╡ 7418aa1c-cc4d-47cc-a132-26f01c0aa761
md"- Selected event = $evt0 (use window to change selection)"

# ╔═╡ 4c6aa884-e7f2-4923-84f9-e78869670e1e
evt = jp.JPetalo.select_event(pdf.total_charge, evt0)

# ╔═╡ 9e48d056-0f4e-4e5b-9a8d-3d0e365a798a
sids = evt[!,:sensor_id]

# ╔═╡ ac738ea2-d2fd-4f0b-9b97-dc1745cb8e22
md"- Function `sipm_xyzq(evt, sensor_xyz)` takes a DataFrame represeting the event and the database with the position of the SiPMs (also a DataFrame) to return a hit-DataFrame --a hit is a structure (x,y,z,q)--"

# ╔═╡ 6feca4b1-9af3-481c-9559-74b604513b07
hitdf = jp.JPetalo.sipm_xyzq(evt, pdf.sensor_xyz)

# ╔═╡ df6aaf5c-5cd8-4d2f-a4a1-f367c4164c5b
length(hitdf.q)

# ╔═╡ aa9584a3-9d41-493a-b4a6-0862e3354f80
evt.charge

# ╔═╡ db157d7e-70e5-4d20-8fab-92454b5c2e09
@test hitdf.q == evt.charge

# ╔═╡ c471cec7-a1cb-4ffe-ac01-8eaf11e2ed1e
md"`plot_xyzq(hdf)` takes a hit DataFrame and returns the scatter plot of the SiPM hits and the tan(y/x) which allows to discriminate one gamma from the other "

# ╔═╡ 7019b621-6ef5-46e8-abb2-f748f005201d
pxyq, pxy, pxz, pyz = jp.JPetalo.plot_xyzq(hitdf, 100.);

# ╔═╡ 5a75154c-bde7-4a2f-af17-63b524ad5958
plot(pxyq,pxy,pxz,pyz, layout = (2, 2), aspect_ratio=:equal,size = (1400, 1000), legend = false,  fmt = :png)


# ╔═╡ 00c704f6-c518-4a05-b49b-b78972231ffe
md"- Histogramming the energy of the SiPMs one can see a large peak at 1 pes. This is the background (diffused light), and suggests a cut at 1 pes"

# ╔═╡ 90777a77-4523-4496-a21f-42338f7e3079
histogram(hitdf.q[hitdf.q.<10], bins=30)

# ╔═╡ 0d8a7dd5-6529-4d0d-a61b-31893cf92262
@bind ecut NumberField(1.0:100.0; default=5.0)

# ╔═╡ b95a86e9-3503-47b3-973a-ee168616b695
typeof(ecut)

# ╔═╡ 5a7d478a-9a79-415d-b99f-d0c548cebeb7
md" - ecut =$ecut pes (can be changed using window)"

# ╔═╡ 4be9f445-e372-4999-a4cb-67565559b6b5
evtQ1 = evt[evt.charge.>ecut,:]

# ╔═╡ 493721a5-7875-48a0-bc46-735021caa292
hitQdf = jp.JPetalo.sipm_xyzq(evtQ1, pdf.sensor_xyz)

# ╔═╡ 5500660e-d169-4fa4-b53f-c32dfdb0e15a
histogram(hitQdf.q, bins=30)

# ╔═╡ 0e880066-fc28-4fde-97e2-3b4e64f1d5b3
size(hitQdf)

# ╔═╡ 0c74c3cc-bcc0-4da8-ad7d-fc00e4cb5c20
function ksipmsel(hdf::DataFrame, ka)
	return hdf[(ka.==2), :], hdf[(ka.==1), :]
end

# ╔═╡ b7625433-4509-422b-9cdd-83b9984f43da
md"- Repeating the plots ater cut at $ecut pes shows a much cleaner distribution"

# ╔═╡ 8e882ad9-2a98-438b-883f-6141c67375cb
pxyqQ1,pxyQ1,pxzQ1,pyzQ1 = jp.JPetalo.plot_xyzq(hitQdf, 100.);

# ╔═╡ dd577b58-33ff-4852-81c2-1d53918234f9
plot(pxyqQ1,pxyQ1,pxzQ1,pyzQ1, layout = (2, 2), aspect_ratio=:equal,size = (1400, 1000), legend = false,  fmt = :png)

# ╔═╡ 35132119-c793-4bfe-b228-7a017ce7789d
md"- We can now select the hits with positive and negative phi, which define the clusters of the gammas"

# ╔═╡ 2af988f3-4754-4de8-833b-a5bc57f0381d
hqpdf, hqndf = jp.JPetalo.sipmsel(hitQdf);

# ╔═╡ da5306f6-9bea-4cc3-9bf9-767b0908fb69
md"- And compute the baricenters, which define the gamma position"

# ╔═╡ 7f277fff-efaa-4502-8a58-45c4d4a514f5
bp = jp.JPetalo.baricenter(hqpdf)

# ╔═╡ dbd61080-79d6-4e64-8a87-a1ed243ff503
bn = jp.JPetalo.baricenter(hqndf)

# ╔═╡ 75c15d2b-c314-493a-a773-cdffc61b43f2
C=transpose([bp.x bp.y bp.z; bn.x bn.y bn.z])

# ╔═╡ 72b3c6a8-cd90-40de-b1b0-2907588ecc92
md"- Finally we can draw the three proyections of the baricenter, together with the LORs that connect them"

# ╔═╡ 9a07dac6-f038-429c-a51a-a4237532fe82
sxy, syz, sxz = jp.JPetalo.plot_barycenter(bp,bn, pxyQ1,pxzQ1,pyzQ1, 0.1);

# ╔═╡ e7eda113-32ad-47b9-8a74-10f759165e16
plot(sxy,syz,sxz, layout = (1, 3), aspect_ratio=:equal,size = (1400, 1000),legend = false,  fmt = :png)

# ╔═╡ 7941d91e-9807-48b7-b454-f8f192a2695c
md"## Reconstruction of LORS"

# ╔═╡ d9aeb478-23a5-4b5d-a993-873223cbd224
typeof(pdf)

# ╔═╡ d0193465-d2f0-4232-994d-4f3a435ed42b
typeof(th1)

# ╔═╡ 081194bf-31a5-47fb-8fe2-79c2e54ae88b
typeof(ecut)

# ╔═╡ 7cfe5f8e-c951-49d0-89d7-7211c24c2f21
typeof(th1.event_id)

# ╔═╡ f8e5dd82-0a49-4d84-8a29-3cb221fce0a5
QM, BP, BN, rBP, rBN = jp.JPetalo.reco_lor(pdf, th1.event_id, ecut)

# ╔═╡ 83ebaade-d34f-4631-a2a3-3e1012113d59


# ╔═╡ 54743f18-0f96-415f-a821-9743b33953be
histogram(QM, bins=200, xlim=[10.,1500.])

# ╔═╡ ae9283ae-fa33-4057-beb8-b49416ea90bb
minimum(QM)

# ╔═╡ 521b5044-6492-468b-8e49-57fd196e40d9
BP

# ╔═╡ 7e933768-d917-4daa-9857-c4b1e72fd0c4
size(BP)

# ╔═╡ 83ed2f0e-3e88-42fe-91f6-a51a21881d8b
rBP

# ╔═╡ ec1aab0e-fb78-42af-b1f4-e41ce9e19389
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

# ╔═╡ a2c0da22-3df5-4b92-9f1d-562c4b5c1087
bsxy, bsxz, bsyz = jp.JPetalo.plot_lors_barycenter(BP, BN, 1500.) 

# ╔═╡ 60e090f2-b35d-4333-9464-6a9d9ac0084d
rsxy, rsxz, rsyz = jp.JPetalo.plot_lors_barycenter(rBP, rBN, 1500.) 

# ╔═╡ 8e51931c-cd9e-4056-816e-9a8e7a1f59a8
plot(bsxy, rsxy, layout= (1, 2), aspect_ratio=:equal,size = (1400, 1000),legend = false,  fmt = :png)

# ╔═╡ edc1f1e2-9d55-4c00-a2af-08571c7f38ab
plot(bsxz, rsxz, layout = (2, 1),  aspect_ratio=:equal, size = (1600, 800), ylim=(-100.,100.), legend = false,  fmt = :png)

# ╔═╡ 382f01d9-6041-47ac-8ed6-b6cd1ea092c4
plot(bsyz, rsyz, layout = (2, 1),  aspect_ratio=:equal, size = (1600, 800), ylim=(-100.,100.), legend = false,  fmt = :png)

# ╔═╡ 71acb562-bcfd-40e7-a3b9-dcffe00cd557
plot(bsxy, rsxy, layout= (1, 2), aspect_ratio=:equal,size = (1400, 1000),legend = false,  fmt = :png)

# ╔═╡ 91ed4b00-54a4-4345-8212-ed25c27fe3aa
plot(ptsxy, bsxy, layout= (1, 2), aspect_ratio=:equal,size = (1400, 1000),legend = false,  fmt = :png)

# ╔═╡ 4ce431ab-aac3-40c5-a15a-890d46d4b501
plot(ptsxz, bsxz, layout= (2, 1), aspect_ratio=:equal, size = (1600, 800), ylim=(-100.,100.),legend = false,  fmt = :png)

# ╔═╡ 256eb45e-88b2-41d3-8209-41e770bb9a11
plot(ptsyz, bsyz, layout= (2, 1), aspect_ratio=:equal, size = (1600, 800), ylim=(-100.,100.),legend = false,  fmt = :png)

# ╔═╡ 1709ae5d-b256-4ab6-8c16-5edaa354f867
md"## Types"

# ╔═╡ ccf11f8d-5bcb-40e1-8b6d-d88812c41f88
md"""
	struct Hit
	x::Float64
	y::Float64
	z::Float64
	q::Float64
	

Represent a (high level) SiPM hit
"""

# ╔═╡ 35d0f27d-0f97-401a-9a82-9e176cf62fa2
struct Hit
	x::Float64
	y::Float64
	z::Float64
	q::Float64
end

# ╔═╡ 56b34074-ac36-476b-b3d9-9057fad68693
md"# Functions"

# ╔═╡ a124b244-f873-4847-b5de-a18b80975e80
md"""
	radius(x::Number, y::Number) = sqrt(x^2 + y^2)
"""

# ╔═╡ b7a9953e-4392-4762-a6f4-979d47426639
begin
	radius(x::Number, y::Number) = sqrt(x^2 + y^2)
	radius(x::Float64, y::Float64) = sqrt(x^2 + y^2)
end

# ╔═╡ 567ca9fc-ce60-4b30-9246-07d999cd7654
md"""
	find_max_xy(df, xc, yc)

Return ymax and x such that ymax = f(x).

**Description:**

In a DataFrame one has often "XY" variables, that is, a pair of columns
"X" and "Y" which represent correlated variables (e.g intensity and wavelength).
In such cases one often wants the XY maximum, that is, finding the maximum
of "Y" (call it ymax) and the corresponding value in "X"
(x for which y is ymax). This corresponds, for example, to the wavelength at
which the intensity is maximal.


**Arguments:**
- `df::DataFrame`: data frame holding the data.
- `xc::String`: the X column.
- `yc::String`: the Y column.
"""

# ╔═╡ 12e87e3c-5521-4d94-b7f6-009abf5e96f1
function find_max_xy(df::DataFrame, xc::String, yc::String)
	ymax, imax = findmax(df[!, yc])
	x_ymax = df[imax, xc]
	return ymax, x_ymax
end

# ╔═╡ 6e112f16-c32c-491a-9b68-180a57aace8f
find_max_xy(evt,"sensor_id", "charge")

# ╔═╡ c91baad0-7c07-4b40-a99e-3c773b18b065
scatter(evt.sensor_id, evt.charge, marker=:circle, leg=false)

# ╔═╡ 982329e8-61df-4bb9-b5fa-82f5c50edafe
function test_find_max_xy(evt)
	qmax, iqmax = findmax(evt.charge)
	simax = evt.sensor_id[iqmax]
	qm, sim = find_max_xy(evt,"sensor_id", "charge")
	return qm ≈ qmax && sim ≈ simax
end

# ╔═╡ 3f26bf3b-a006-4956-9674-09b157c247c9
@test test_find_max_xy(evt)

# ╔═╡ 300b1614-741d-4720-b2bc-01539423a163
md"""
	find_xyz_sipm_qmax(hitdf::DataFrame)
Return the coordinates of the SiPM with maximum charge
"""

# ╔═╡ 3872788b-62a9-4733-8d88-b4af0e59c423
function find_xyz_sipm_qmax(hitdf::DataFrame)
	qxmax, xmax = find_max_xy(hitdf,"x", "q")
	qymax, ymax = find_max_xy(hitdf,"y", "q")
	qzmax, zmax = find_max_xy(hitdf,"z", "q")
	return Hit(xmax,ymax,zmax,qxmax)
end

# ╔═╡ d09f2c56-4ac1-47b4-b5d2-534a4c7b5301
function test_find_xyz_sipm_qmax(hitdf)
	qxmax, xmax = find_max_xy(hitdf,"x", "q")
	qymax, ymax = find_max_xy(hitdf,"y", "q")
	qzmax, zmax = find_max_xy(hitdf,"z", "q")
	@test qxmax ≈ qymax
	@test qxmax ≈ qzmax
end

# ╔═╡ 44af890b-5d5b-4e57-b284-0a7729674466
test_find_xyz_sipm_qmax(hitdf)

# ╔═╡ 1cdcefbf-99e1-43b3-9587-b9b5bfee40c4
md" - `simax` is the Hit with the coordinates and charge of the SiPM of max charge"

# ╔═╡ e694f819-3099-4a90-9b44-e2151112941c
simax = find_xyz_sipm_qmax(hitdf::DataFrame)

# ╔═╡ 0a2569da-b161-428b-b092-722e01255971
md"""
	xyz_dot(hitdf::DataFrame, simax::Hit)
Return the dot product between each SiPM in the event and the SiPM of max charge
"""

# ╔═╡ 735c6816-a7fc-457f-ae8b-95886fcc84bb
function xyz_dot(hitdf::DataFrame, simax::Hit)
	xyzmax = [simax.x, 	simax.y, simax.z]
	xyzm_dot = dot(xyzmax, xyzmax)
	return [dot(Array(hitdf[i,1:3]), xyzmax) /xyzm_dot  for i in 1:nrow(hitdf)]
end

# ╔═╡ 81685b61-d54e-4c6b-a087-c53a7d0f62e3
md"Plot `xyz_dot` showing that the SiPM are distributed around zero and between -π/2 and π/2"

# ╔═╡ 68847964-7edf-4dba-9cd0-2e9197fdf5a9
npr = xyz_dot(hitdf, simax)

# ╔═╡ bbe9358e-1d79-4bd2-b003-d49831e21cfc
histogram(npr)

# ╔═╡ 072d42e4-9924-4d2d-8302-c2cd34ec9ee3
@test maximum(npr) < π/2

# ╔═╡ 2c52dc80-04ca-48e3-acfc-3d6f63acb323
@test minimum(npr) > -π/2 

# ╔═╡ 8d38fc77-7801-42b6-91d4-d56e3f657f1b
md"""
	sipmsel(hdf::DataFrame)
Return two data frames, separating the SiPMs in the phi angle relative to the SiPM of max charge. 
"""

# ╔═╡ c9111054-494c-4cc7-96cd-26d4c522da1b
function sipmsel(hdf::DataFrame)
	simax = find_xyz_sipm_qmax(hdf)
	npr   = xyz_dot(hdf, simax)
	mask =[n>0 ? true : false for n in npr]
	return hdf[(npr.>0), :], hdf[(npr.<0), :]
end

# ╔═╡ 050e991d-707f-43c2-80cd-d8ce570ac360
function fphi(hdf::DataFrame)
	return atan.(hdf.y,hdf.x)
end

# ╔═╡ 33192e19-2a20-4a43-ba36-2a5d40d009cf
hitp, hitn  = sipmsel(hitdf)

# ╔═╡ 12044636-6967-424f-977f-c11e8df0b2e1
scatter(fphi(hitdf), leg=false)

# ╔═╡ 7f8dd559-5c77-4924-9eb9-dfb6eabbed09
md"""
	select_by_column_value(df::DataFrame, column::String, value)

Select the elements of a dataframe in terms of value.
"""

# ╔═╡ 91e8eb43-3dde-4472-9cdf-07f86fde14a6
function select_by_column_value(df::DataFrame, column::String, value)
	mask = df[!,column].==value
	return df[mask,:]
end

# ╔═╡ 459e546d-825a-431e-a676-874a9bb926c4
md"""
	select_by_index(df::DataFrame, column::String, value::Int64)

Select the elements of a dataframe in terms of an index
"""

# ╔═╡ c925cb36-358b-4e4a-82a8-11f06d162a11
select_by_index(df::DataFrame, column::String, value::Int64) = select_by_column_value(df::DataFrame, column::String, value)

# ╔═╡ 3ae4d41f-624a-4f24-8d1f-d531ada29407
md"- We can test that the positions of the sipms corresponding to a given index are correctly obtained"

# ╔═╡ 0b7a4788-9979-465b-bf23-b25f9ffc89f6
select_by_index(pdb, "SensorID", 1000)

# ╔═╡ daf2a651-519d-4861-865c-ca696ec4f419
md"""
	sipm_pos(dbdf::DataFrame, index::Int64)
Take the sipm database and the index of a sipm and returns its position
"""

# ╔═╡ 0e32b95a-587f-4f30-9f47-49db52fedf9e
function sipm_pos(dbdf::DataFrame, index::Int64)
	return Array(select_by_index(dbdf, "SensorID", index)[1,2:end])
end

# ╔═╡ ba09133d-9fd7-43b9-b18b-2be8e95bc56e
@test sipm_pos(pdb, 1000) ≈ Array(select_by_index(pdb, "SensorID", 1000)[1,2:end])

# ╔═╡ 16c3f02f-7ce8-4664-881b-9a93dfd18f54
md"- Example: Function applied over an index"

# ╔═╡ 4ca52cab-b70e-4f56-970a-be78fef81a27
x1,y1, z1 = sipm_pos(pdb, 77194)

# ╔═╡ 708c295b-5fb3-4e45-8674-389356576139
md"- Function applied over all indexes"

# ╔═╡ 86a68e86-8ad1-4c09-8bd7-178106d35cf3
begin
	XYZ = sipm_pos.((pdb,), pdb.SensorID)
	x = [XYZ[i][1] for i in 1:length(XYZ)]
	y = [XYZ[i][2] for i in 1:length(XYZ)]
	z = [XYZ[i][3] for i in 1:length(XYZ)]
	r = radius.(x, y)
end

# ╔═╡ 5c4f390d-7224-49db-b16d-2fd3339bd47b
function get_hits_as_matrix(hitdf)
	f = @formula(0 ~ x + y + z)
	f2 = apply_schema(f, schema(f, hitdf))
	resp, pred = modelcols(f2, hitdf)
	return transpose(pred)
end

# ╔═╡ 1fe2246a-4944-4638-9c0a-309eec3c6452
Mhits = get_hits_as_matrix(hitQdf)

# ╔═╡ 4ac855d4-396b-4161-ab36-0b6fae24c98d
kr = kmeans(Mhits, 2)

# ╔═╡ 99fd2114-a8f1-45f0-aa52-bc715990edd2
begin
	ka = assignments(kr) # get the assignments of points to clusters
	kc = counts(kr) # get the cluster sizes
	kM = kr.centers # get the cluster centers
end

# ╔═╡ fe7a6181-c975-435f-8c34-1e7f9e5f85c3
ka

# ╔═╡ ef61b82b-f616-43af-af9e-6f8d6451bd6e
size(ka)

# ╔═╡ 88639225-8441-4f21-851b-9b75a783b1ce
typeof(ka)

# ╔═╡ 6090917b-50f2-452f-9fbd-deea914d3de6
kc

# ╔═╡ 4a0d90d7-3a97-46ce-926e-f906459108ad
kc[1] + kc[2]

# ╔═╡ 5d3e5e32-1d3f-4551-90be-68bf5d78920e
hq2df, hq1df = ksipmsel(hitQdf, ka)

# ╔═╡ 2a4af456-47c3-4330-8f8b-6d3812150e10
b2 = jp.JPetalo.baricenter(hq2df)

# ╔═╡ 23d0f15e-d54f-4db1-ad35-c4de300db9ba
b1 = jp.JPetalo.baricenter(hq1df)

# ╔═╡ 4f1e223d-a51a-4ec4-afca-358a1f0d0139
kM

# ╔═╡ b215f8f2-2a09-467d-b831-50daad8ca542
R = kmeans!(Mhits, C)

# ╔═╡ b73f5a70-a0ce-4ebb-8dfc-b006f6acd242
rk = R.centers

# ╔═╡ fe513e3f-d39c-4e22-af36-f66f353623d2
rbp = jp.JPetalo.Hit(rk[1,1],rk[2,1], rk[3,1], 0.)

# ╔═╡ a41bfeb3-7b1f-4706-922c-7abfa9aadb55
rbn = jp.JPetalo.Hit(rk[1,2],rk[2,2], rk[3,2], 0.)

# ╔═╡ 9280feb6-636b-49b7-b886-28e93deddd2f
function reco_lor(pdf, th1, ecut)
	function ksipmsel(hdf::DataFrame, ka)
		return hdf[(ka.==2), :], hdf[(ka.==1), :]
	end
	
	BP = []
	BN = []
	rBP = []
	rBN = []
	QM = []
	for evt0 in th1.event_id
		evt = jp.JPetalo.select_event(pdf.total_charge, evt0)
		evtQ1 = evt[evt.charge.>ecut,:]
		hitQdf = jp.JPetalo.sipm_xyzq(evtQ1, pdf.sensor_xyz)
		qm = maximum(hitQdf.q)
		
		Mhits = get_hits_as_matrix(hitQdf)
		kr = kmeans(Mhits, 2)
		ka = assignments(kr) # get the assignments of points to clusters
		kc = counts(kr) # get the cluster sizes
		rk = kr.centers # get the cluster centers
		
		hq2df, hq1df = ksipmsel(hitQdf, ka)
		b1 = jp.JPetalo.baricenter(hq1df)
		b2 = jp.JPetalo.baricenter(hq2df)
		
		rb1 = jp.JPetalo.Hit(rk[1,1],rk[2,1], rk[3,1], b1.q)
		rb2 = jp.JPetalo.Hit(rk[1,2],rk[2,2], rk[3,2], b2.q)
		
		push!(BP, b1)
		push!(BN, b2)
		push!(rBP, rb1)
		push!(rBN, rb2)
		push!(QM, qm)

	end
	return QM, BP, BN, rBP, rBN
end

# ╔═╡ 08746abc-529b-4bef-baf3-24a01f215821
scatter(x, y)

# ╔═╡ 64fe61ae-0909-4aa5-8f6c-223faeb73d34
histogram(z)

# ╔═╡ 3b051028-9e0f-4c78-99cd-73d6ed0e61b7
@test mean(r) ≈ 409.6000094986399

# ╔═╡ 389b5b10-0bad-483e-b02a-3dc15e231c29
mean(r)

# ╔═╡ 6d992b8f-210a-468c-9a72-5bda8564f619
@test mean(z) ≈0

# ╔═╡ 138db3d0-4920-4405-949f-d1ff3308c696
md"""
	select_event(dbdf::DataFrame, index::Int64)

Take the event dataframe and the index of an event and returns a data frame which selects that particular event

"""

# ╔═╡ dfcfd7e7-b6f8-423c-a177-be6360a031c1
function select_event(dbdf::DataFrame, index::Int64)
	return select_by_index(dbdf, "event_id", index)[:,2:end]
end

# ╔═╡ 86078e18-ccbf-40cd-9b4c-08f1a5469147
md"""
	sipm_xyzq(evt::DataFrame, pdb::DataFrame)

Take an event dataframe and a dataframe with the sipm positions and returns a hit-dataframe, with (x,y,z,q) for each sipm
"""

# ╔═╡ dc324771-35b0-4410-a0d5-8582d8975d3a
function sipm_xyzq(evt::DataFrame, pdb::DataFrame)
	sids = evt[!,:sensor_id]
	pos = sipm_pos.((pdb,),sids)
	x = [p[1] for p in pos]
	y = [p[2] for p in pos]
	z = [p[3] for p in pos]
	q = evt[!,:charge]
	return DataFrame(x=x,y=y,z=z,q=q)
end

# ╔═╡ 7623540b-7892-45b0-bd28-39c2fedf620d
md"""
	baricenter(hdf::DataFrame)

Given a Hit DataFrame (a dataframe with columns x,y,z,q) return its barycenter
"""

# ╔═╡ f743a6f0-7f28-4bb7-be72-b25424879312
function baricenter(hdf::DataFrame)
	function xq(hdf::DataFrame, pos::String)
		return sum(hdf[!,pos] .* hdf.q) / qt
	end
	qt = sum(hdf.q)
	return Hit(xq(hdf, "x"), xq(hdf, "y"), xq(hdf, "z"), qt)
end

# ╔═╡ d859296b-b081-4e3b-aa9c-08c8f6fc99cd
md"## Plot functions"

# ╔═╡ 02434a19-c1cf-4c60-9457-928d64e3419b
md"""
	plot_xyzq(hdf::DataFrame)
Take a hit dataframe and produce 4 scatter plots showing the clustering of the SiPMs and the phi distribution showing how the clusters distribute at opposite phi values 
"""

# ╔═╡ 625b20ac-1c95-4206-85bb-14705a4028f4
function plot_xyzq(hdf::DataFrame)
	pxyq = scatter(hdf.x,hdf.y,hdf.q, legend=false)
	pxy = scatter(hdf.x,hdf.y, marker_z = hdf.q, markersize = hdf.q/3,  color = :jet,
		         legend=false)
	xlabel!("x")
	ylabel!("y")
	pxz = scatter(hdf.x,hdf.z, marker_z = hdf.q, markersize = hdf.q/3,  color = :jet,
		         legend=false)
	xlabel!("x")
	ylabel!("z")
	pyz = scatter(hdf.y,hdf.z, marker_z = hdf.q, markersize = hdf.q/3,  color = :jet,
		         legend=false)
	xlabel!("y")
	ylabel!("z")
	pphi = scatter(fphi(hdf), leg=false)
	pphi = plot(atan.(hdf.y,hdf.x), shape  = :circle, legend=false)
	#xlabel!("SiPM number")
	ylabel!("tan(y/x)")
	hphi = histogram(atan.(hdf.y,hdf.x), bins=20)
	xlabel!("tan(y/x)")
	return pxyq,pxy,pxz,pyz, pphi, hphi
end

# ╔═╡ afac8cb0-83d3-4be8-90a4-329b15d1b614
md"""
	plot_barycenter(bp,bn, pxy, pxz, pyz)
Plots the barycenter of the SiPm clusters together with the LOR that connects them
"""

# ╔═╡ 856e70b5-0f27-45c6-afa2-b86cfd19f1bb
function plot_barycenter(bp,bn, pxy, pxz, pyz)
	lsxy = LineSegment([bp.x,bp.y],[bn.x,bn.y])
	sxy = scatter!(pxy,[bp.x,bn.x],[bp.y,bn.y], 
		marker_z = [bp.q,bn.q], markersize = 5,  color = :jet, legend=false)
	sxy = plot!(sxy,lsxy)
	syz = scatter!(pyz,[bp.y,bn.y],[bp.z,bn.z], 
		marker_z = [bp.q,bn.q], markersize = 5, color = :jet,legend=false)
	lsyz = LineSegment([bp.y,bp.z],[bn.y,bn.z])
	syz  = plot!(syz,lsyz)
	sxz = scatter!(pxz,[bp.x,bn.x],[bp.z,bn.z], 
		marker_z = [bp.q,bn.q], markersize = 5, color = :jet,legend=false)
	lsxz = LineSegment([bp.x,bp.z],[bn.x,bn.z])
	sxz  = plot!(sxz,lsxz)
	return sxy, syz, sxz
end

# ╔═╡ Cell order:
# ╠═79cfd2fc-9046-11eb-2b13-1b877d57645d
# ╠═d16a3879-acb4-4cce-9b87-faadf4abfece
# ╠═4ee98da9-f3ce-4782-9538-f878f27ed9f7
# ╠═2e5dd476-4e23-427e-8e72-a6af044eb397
# ╠═26a0d201-f36c-4bcf-9219-a2f22c500598
# ╠═3207f446-8643-11eb-37ba-c9aec47fcb8f
# ╠═5115917a-8644-11eb-19fc-0528741ca75d
# ╠═fc8b79a2-8728-11eb-2da7-e3ffa3ceef08
# ╠═ac019a72-dca3-4e7e-9ffe-3cfbede4106e
# ╠═aa9f1b97-21cb-4592-8f47-8549f1305da2
# ╠═0068d015-f292-4bf7-82f9-c6c0115f96e2
# ╠═68e738e7-88bd-41c2-89e4-594f07d64ddc
# ╠═0f2f4c78-8729-11eb-2bab-27812ce8c47e
# ╠═5ee27d52-86fd-11eb-365e-9f2b1a095575
# ╠═621ec96c-86fd-11eb-1c41-379cc17180dc
# ╠═9b853f27-4288-42a6-8f12-ca004e1773b7
# ╠═26225a47-c7aa-4319-89f2-1abbc5c1116a
# ╠═cf89b973-7b0f-483d-8ce9-ba426f1df2a6
# ╠═0b1cedc7-8ada-45a5-adef-fbae794dee3e
# ╠═5fb08873-1ca6-44f2-b68a-438fea6007ed
# ╠═383048a4-12e4-492f-bcbd-a83c6e7dc7aa
# ╟─80542fd1-843e-4d78-9bd3-169c1d6a9672
# ╠═edbcd276-6e0d-41f3-91d8-131b0fc7486b
# ╠═fdff3aeb-834c-4726-8e6c-a35e6b9d46d1
# ╠═bc16f5c1-863e-47de-87aa-838578157188
# ╠═5b8468a8-4465-4708-afd0-9c7495b3a8a3
# ╠═8db332bb-9a40-4bdf-b018-64f6a35f6bdf
# ╠═b8aacd81-9d8f-4105-a214-220be62b3e43
# ╠═42e1caa1-4aba-44e7-ab32-4538662bcf29
# ╠═fc41144a-ba14-4543-9523-c5450744e126
# ╠═a2f23d32-150b-4bb7-96d0-79f687747e40
# ╠═6fd444d3-03f3-4b56-a02d-cddba85ca96f
# ╠═06208989-8940-4cd5-8f59-c5e7c167b97b
# ╠═5c4b2813-991e-4908-833c-195309dd95f1
# ╠═bfc9e52f-7889-48fc-a6a8-55b8844420e0
# ╠═abfb7e89-99dd-4e2a-8b5f-ec8fefb8fdc3
# ╟─766c0b57-9a16-44d1-ab5f-e7c1d39a5176
# ╠═48a21c9c-cfc3-4d23-b1d9-cb9e47753bb4
# ╠═361424ff-0d5d-4150-8928-928a184db66f
# ╠═2231ac80-f3cf-446b-b6dd-ba02efbd29b7
# ╠═9e792d68-79a7-4418-b96f-13c32be0e20c
# ╠═868a8593-f112-43ed-af79-88d972fd2791
# ╠═3e514c2f-03b9-42da-8317-0c80309066a9
# ╠═e025a1d1-d8ef-4ee7-8c7a-a58eb4d6eee2
# ╠═cbd4bcba-1429-4bcf-8fec-20316b0f18e2
# ╠═b8d08354-9346-4e70-ab23-9a7927524099
# ╠═0a7cfeff-e1e4-4662-bc59-5088e95749b6
# ╠═6c355936-7d4a-456c-bddb-0de22d3144fe
# ╠═bb427097-45b5-42a3-9c4b-ab48913187bd
# ╠═ce932b08-2e26-4527-94f1-da3c49a5433b
# ╠═99500b8e-3559-495e-83e9-42dde08173b6
# ╠═9963d767-a54b-40db-ba6a-a3d64f6dcbd7
# ╠═df6d7369-88ae-4536-b66c-369d56291bce
# ╠═21d2c293-ec33-48a7-83f4-8cd6df8c0a70
# ╠═4d24d7fa-28f7-4d0a-9889-f0790ae6fe53
# ╠═302b25f2-cd53-4918-983d-dc37399507a3
# ╠═178d5a6c-2fb8-49cc-884b-414f86ceaad6
# ╠═2694fc3f-2883-4af7-ba40-6a5b41ce1147
# ╠═f7b06fae-650c-4ea2-9a12-4a28b775f845
# ╠═f1d2767f-c729-4f22-8e83-5cd3161f7d93
# ╠═8b050930-1c99-4613-aa95-90f40743acab
# ╠═f4f30742-121d-49ea-ad14-df9650bf6c6b
# ╠═801ded29-0d7f-4407-94e8-7fffb95613b8
# ╠═8a6c603d-a33b-440b-b5bd-9b7bef37b2b5
# ╠═f3b71ecf-fa44-480d-bb7d-612edc6ff15e
# ╟─7418aa1c-cc4d-47cc-a132-26f01c0aa761
# ╠═4c6aa884-e7f2-4923-84f9-e78869670e1e
# ╠═9e48d056-0f4e-4e5b-9a8d-3d0e365a798a
# ╠═ac738ea2-d2fd-4f0b-9b97-dc1745cb8e22
# ╠═6feca4b1-9af3-481c-9559-74b604513b07
# ╠═df6aaf5c-5cd8-4d2f-a4a1-f367c4164c5b
# ╠═aa9584a3-9d41-493a-b4a6-0862e3354f80
# ╠═db157d7e-70e5-4d20-8fab-92454b5c2e09
# ╟─c471cec7-a1cb-4ffe-ac01-8eaf11e2ed1e
# ╠═7019b621-6ef5-46e8-abb2-f748f005201d
# ╠═5a75154c-bde7-4a2f-af17-63b524ad5958
# ╟─00c704f6-c518-4a05-b49b-b78972231ffe
# ╠═90777a77-4523-4496-a21f-42338f7e3079
# ╠═0d8a7dd5-6529-4d0d-a61b-31893cf92262
# ╠═b95a86e9-3503-47b3-973a-ee168616b695
# ╟─5a7d478a-9a79-415d-b99f-d0c548cebeb7
# ╠═4be9f445-e372-4999-a4cb-67565559b6b5
# ╠═493721a5-7875-48a0-bc46-735021caa292
# ╠═5500660e-d169-4fa4-b53f-c32dfdb0e15a
# ╠═5c4f390d-7224-49db-b16d-2fd3339bd47b
# ╠═1fe2246a-4944-4638-9c0a-309eec3c6452
# ╠═4ac855d4-396b-4161-ab36-0b6fae24c98d
# ╠═99fd2114-a8f1-45f0-aa52-bc715990edd2
# ╠═fe7a6181-c975-435f-8c34-1e7f9e5f85c3
# ╠═ef61b82b-f616-43af-af9e-6f8d6451bd6e
# ╠═88639225-8441-4f21-851b-9b75a783b1ce
# ╠═0e880066-fc28-4fde-97e2-3b4e64f1d5b3
# ╠═6090917b-50f2-452f-9fbd-deea914d3de6
# ╠═4a0d90d7-3a97-46ce-926e-f906459108ad
# ╠═0c74c3cc-bcc0-4da8-ad7d-fc00e4cb5c20
# ╠═5d3e5e32-1d3f-4551-90be-68bf5d78920e
# ╟─b7625433-4509-422b-9cdd-83b9984f43da
# ╠═8e882ad9-2a98-438b-883f-6141c67375cb
# ╠═dd577b58-33ff-4852-81c2-1d53918234f9
# ╟─35132119-c793-4bfe-b228-7a017ce7789d
# ╠═2af988f3-4754-4de8-833b-a5bc57f0381d
# ╟─da5306f6-9bea-4cc3-9bf9-767b0908fb69
# ╠═7f277fff-efaa-4502-8a58-45c4d4a514f5
# ╠═dbd61080-79d6-4e64-8a87-a1ed243ff503
# ╠═2a4af456-47c3-4330-8f8b-6d3812150e10
# ╠═23d0f15e-d54f-4db1-ad35-c4de300db9ba
# ╠═4f1e223d-a51a-4ec4-afca-358a1f0d0139
# ╠═75c15d2b-c314-493a-a773-cdffc61b43f2
# ╠═b215f8f2-2a09-467d-b831-50daad8ca542
# ╠═b73f5a70-a0ce-4ebb-8dfc-b006f6acd242
# ╠═fe513e3f-d39c-4e22-af36-f66f353623d2
# ╠═a41bfeb3-7b1f-4706-922c-7abfa9aadb55
# ╟─72b3c6a8-cd90-40de-b1b0-2907588ecc92
# ╠═9a07dac6-f038-429c-a51a-a4237532fe82
# ╠═e7eda113-32ad-47b9-8a74-10f759165e16
# ╠═7941d91e-9807-48b7-b454-f8f192a2695c
# ╠═d9aeb478-23a5-4b5d-a993-873223cbd224
# ╠═d0193465-d2f0-4232-994d-4f3a435ed42b
# ╠═081194bf-31a5-47fb-8fe2-79c2e54ae88b
# ╠═9280feb6-636b-49b7-b886-28e93deddd2f
# ╠═7cfe5f8e-c951-49d0-89d7-7211c24c2f21
# ╠═f8e5dd82-0a49-4d84-8a29-3cb221fce0a5
# ╠═83ebaade-d34f-4631-a2a3-3e1012113d59
# ╠═54743f18-0f96-415f-a821-9743b33953be
# ╠═ae9283ae-fa33-4057-beb8-b49416ea90bb
# ╠═521b5044-6492-468b-8e49-57fd196e40d9
# ╠═7e933768-d917-4daa-9857-c4b1e72fd0c4
# ╠═83ed2f0e-3e88-42fe-91f6-a51a21881d8b
# ╠═ec1aab0e-fb78-42af-b1f4-e41ce9e19389
# ╠═a2c0da22-3df5-4b92-9f1d-562c4b5c1087
# ╠═60e090f2-b35d-4333-9464-6a9d9ac0084d
# ╠═8e51931c-cd9e-4056-816e-9a8e7a1f59a8
# ╠═edc1f1e2-9d55-4c00-a2af-08571c7f38ab
# ╠═382f01d9-6041-47ac-8ed6-b6cd1ea092c4
# ╠═71acb562-bcfd-40e7-a3b9-dcffe00cd557
# ╠═91ed4b00-54a4-4345-8212-ed25c27fe3aa
# ╠═4ce431ab-aac3-40c5-a15a-890d46d4b501
# ╠═256eb45e-88b2-41d3-8209-41e770bb9a11
# ╠═1709ae5d-b256-4ab6-8c16-5edaa354f867
# ╟─ccf11f8d-5bcb-40e1-8b6d-d88812c41f88
# ╠═35d0f27d-0f97-401a-9a82-9e176cf62fa2
# ╟─56b34074-ac36-476b-b3d9-9057fad68693
# ╟─a124b244-f873-4847-b5de-a18b80975e80
# ╠═b7a9953e-4392-4762-a6f4-979d47426639
# ╠═567ca9fc-ce60-4b30-9246-07d999cd7654
# ╠═12e87e3c-5521-4d94-b7f6-009abf5e96f1
# ╠═6e112f16-c32c-491a-9b68-180a57aace8f
# ╠═c91baad0-7c07-4b40-a99e-3c773b18b065
# ╠═982329e8-61df-4bb9-b5fa-82f5c50edafe
# ╠═3f26bf3b-a006-4956-9674-09b157c247c9
# ╠═300b1614-741d-4720-b2bc-01539423a163
# ╠═3872788b-62a9-4733-8d88-b4af0e59c423
# ╠═d09f2c56-4ac1-47b4-b5d2-534a4c7b5301
# ╠═44af890b-5d5b-4e57-b284-0a7729674466
# ╟─1cdcefbf-99e1-43b3-9587-b9b5bfee40c4
# ╠═e694f819-3099-4a90-9b44-e2151112941c
# ╠═0a2569da-b161-428b-b092-722e01255971
# ╠═735c6816-a7fc-457f-ae8b-95886fcc84bb
# ╟─81685b61-d54e-4c6b-a087-c53a7d0f62e3
# ╠═68847964-7edf-4dba-9cd0-2e9197fdf5a9
# ╠═bbe9358e-1d79-4bd2-b003-d49831e21cfc
# ╠═072d42e4-9924-4d2d-8302-c2cd34ec9ee3
# ╠═2c52dc80-04ca-48e3-acfc-3d6f63acb323
# ╠═8d38fc77-7801-42b6-91d4-d56e3f657f1b
# ╠═c9111054-494c-4cc7-96cd-26d4c522da1b
# ╠═050e991d-707f-43c2-80cd-d8ce570ac360
# ╠═33192e19-2a20-4a43-ba36-2a5d40d009cf
# ╠═12044636-6967-424f-977f-c11e8df0b2e1
# ╟─7f8dd559-5c77-4924-9eb9-dfb6eabbed09
# ╠═91e8eb43-3dde-4472-9cdf-07f86fde14a6
# ╟─459e546d-825a-431e-a676-874a9bb926c4
# ╠═c925cb36-358b-4e4a-82a8-11f06d162a11
# ╟─3ae4d41f-624a-4f24-8d1f-d531ada29407
# ╠═0b7a4788-9979-465b-bf23-b25f9ffc89f6
# ╠═ba09133d-9fd7-43b9-b18b-2be8e95bc56e
# ╟─daf2a651-519d-4861-865c-ca696ec4f419
# ╠═0e32b95a-587f-4f30-9f47-49db52fedf9e
# ╠═16c3f02f-7ce8-4664-881b-9a93dfd18f54
# ╟─4ca52cab-b70e-4f56-970a-be78fef81a27
# ╟─708c295b-5fb3-4e45-8674-389356576139
# ╠═86a68e86-8ad1-4c09-8bd7-178106d35cf3
# ╠═08746abc-529b-4bef-baf3-24a01f215821
# ╠═64fe61ae-0909-4aa5-8f6c-223faeb73d34
# ╠═3b051028-9e0f-4c78-99cd-73d6ed0e61b7
# ╠═389b5b10-0bad-483e-b02a-3dc15e231c29
# ╠═6d992b8f-210a-468c-9a72-5bda8564f619
# ╠═138db3d0-4920-4405-949f-d1ff3308c696
# ╠═dfcfd7e7-b6f8-423c-a177-be6360a031c1
# ╟─86078e18-ccbf-40cd-9b4c-08f1a5469147
# ╠═dc324771-35b0-4410-a0d5-8582d8975d3a
# ╟─7623540b-7892-45b0-bd28-39c2fedf620d
# ╠═f743a6f0-7f28-4bb7-be72-b25424879312
# ╟─d859296b-b081-4e3b-aa9c-08c8f6fc99cd
# ╠═02434a19-c1cf-4c60-9457-928d64e3419b
# ╠═625b20ac-1c95-4206-85bb-14705a4028f4
# ╟─afac8cb0-83d3-4be8-90a4-329b15d1b614
# ╠═856e70b5-0f27-45c6-afa2-b86cfd19f1bb
