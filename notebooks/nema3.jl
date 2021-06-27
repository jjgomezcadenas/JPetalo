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
pxyqt,pxyt,pxzt,pyzt = jp.JPetalo.plot_truehits(th1, th2, 101.1);

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

# ╔═╡ 0e880066-fc28-4fde-97e2-3b4e64f1d5b3
size(hitQdf)

# ╔═╡ 6090917b-50f2-452f-9fbd-deea914d3de6
kc

# ╔═╡ 4a0d90d7-3a97-46ce-926e-f906459108ad
kc[1] + kc[2]

# ╔═╡ 0c74c3cc-bcc0-4da8-ad7d-fc00e4cb5c20
function ksipmsel(hdf::DataFrame, ka)
	return hdf[(ka.==2), :], hdf[(ka.==1), :]
end

# ╔═╡ 5d3e5e32-1d3f-4551-90be-68bf5d78920e
hq2df, hq1df = ksipmsel(hitQdf, ka)

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

# ╔═╡ 2a4af456-47c3-4330-8f8b-6d3812150e10
b2 = jp.JPetalo.baricenter(hq2df)

# ╔═╡ 23d0f15e-d54f-4db1-ad35-c4de300db9ba
b1 = jp.JPetalo.baricenter(hq1df)

# ╔═╡ 4f1e223d-a51a-4ec4-afca-358a1f0d0139
kM

# ╔═╡ 75c15d2b-c314-493a-a773-cdffc61b43f2
C=transpose([bp.x bp.y bp.z; bn.x bn.y bn.z])

# ╔═╡ b215f8f2-2a09-467d-b831-50daad8ca542
R = kmeans!(Mhits, C)

# ╔═╡ b73f5a70-a0ce-4ebb-8dfc-b006f6acd242
rk = R.centers

# ╔═╡ fe513e3f-d39c-4e22-af36-f66f353623d2
rbp = jp.JPetalo.Hit(rk[1,1],rk[2,1], rk[3,1], 0.)

# ╔═╡ a41bfeb3-7b1f-4706-922c-7abfa9aadb55
rbn = jp.JPetalo.Hit(rk[1,2],rk[2,2], rk[3,2], 0.)

# ╔═╡ 72b3c6a8-cd90-40de-b1b0-2907588ecc92
md"- Finally we can draw the three proyections of the baricenter, together with the LORs that connect them"

# ╔═╡ 9a07dac6-f038-429c-a51a-a4237532fe82
sxy, syz, sxz = jp.JPetalo.plot_barycenter(bp,bn, pxyQ1,pxzQ1,pyzQ1, 0.1);

# ╔═╡ e7eda113-32ad-47b9-8a74-10f759165e16
plot(sxy,syz,sxz, layout = (1, 3), aspect_ratio=:equal,size = (1400, 1000),legend = false,  fmt = :png)

# ╔═╡ 7941d91e-9807-48b7-b454-f8f192a2695c
md"## Reconstruction of LORS"

# ╔═╡ f8e5dd82-0a49-4d84-8a29-3cb221fce0a5
QM, BP, BN, rBP, rBN = jp.JPetalo.reco_lor(pdf, th1.event_id, ecut)

# ╔═╡ 83ebaade-d34f-4631-a2a3-3e1012113d59


# ╔═╡ 54743f18-0f96-415f-a821-9743b33953be
histogram(QM, bins=200, xlim=[10.,1500.])

# ╔═╡ ae9283ae-fa33-4057-beb8-b49416ea90bb
minimum(QM)

# ╔═╡ a2c0da22-3df5-4b92-9f1d-562c4b5c1087
bsxy, bsxz, bsyz = jp.JPetalo.plot_lors_barycenter(BP, BN, 1500.) ;

# ╔═╡ 60e090f2-b35d-4333-9464-6a9d9ac0084d
rsxy, rsxz, rsyz = jp.JPetalo.plot_lors_barycenter(rBP, rBN, 1500.) ;

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
# ╠═f8e5dd82-0a49-4d84-8a29-3cb221fce0a5
# ╠═83ebaade-d34f-4631-a2a3-3e1012113d59
# ╠═54743f18-0f96-415f-a821-9743b33953be
# ╠═ae9283ae-fa33-4057-beb8-b49416ea90bb
# ╠═a2c0da22-3df5-4b92-9f1d-562c4b5c1087
# ╠═60e090f2-b35d-4333-9464-6a9d9ac0084d
# ╠═8e51931c-cd9e-4056-816e-9a8e7a1f59a8
# ╠═edc1f1e2-9d55-4c00-a2af-08571c7f38ab
# ╠═382f01d9-6041-47ac-8ed6-b6cd1ea092c4
# ╠═71acb562-bcfd-40e7-a3b9-dcffe00cd557
# ╠═91ed4b00-54a4-4345-8212-ed25c27fe3aa
# ╠═4ce431ab-aac3-40c5-a15a-890d46d4b501
# ╠═256eb45e-88b2-41d3-8209-41e770bb9a11
