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

# ╔═╡ 08f6e767-a7cb-4fea-af2d-6345019f6c2b
Pkg.add("Dates")

# ╔═╡ c1c9289b-adf3-4987-8b3c-f1eec3b91388
Pkg.add("Distributions")

# ╔═╡ ad72864c-0d9f-418e-a646-ba231690b31f
using Glob

# ╔═╡ 3207f446-8643-11eb-37ba-c9aec47fcb8f
begin
	using Markdown
	using InteractiveUtils
	using PlutoUI
	using Test
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

# ╔═╡ 0068d015-f292-4bf7-82f9-c6c0115f96e2
using Dates

# ╔═╡ 3ac83772-8efb-4d5f-8639-398d3157e58a
using Plots

# ╔═╡ 37dde9a4-3735-441b-b782-478e364468d6
using Random

# ╔═╡ d215ef1a-817d-4cc7-8c3a-2ba8d5c411cb
using Distributions

# ╔═╡ 79cfd2fc-9046-11eb-2b13-1b877d57645d
md"# NEMA3

- NEMA3 studies
"

# ╔═╡ 4ee98da9-f3ce-4782-9538-f878f27ed9f7
#Pkg.add.(["HDF5"])

# ╔═╡ 2e5dd476-4e23-427e-8e72-a6af044eb397
#Pkg.add.(["Clustering"])

# ╔═╡ 26a0d201-f36c-4bcf-9219-a2f22c500598
#Pkg.add.(["StatsModels"])

# ╔═╡ ba42ec81-5d26-4759-87e1-2bc004a83ca9
#Pkg.add([
			#Pkg.PackageSpec(name="GR", version=v"0.53.0")
		#])

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

# ╔═╡ d29c1b29-4e5e-44aa-8ca0-334e2a74d6d2
n3 = ingredients(srcdir("nema3a.jl"))

# ╔═╡ 0b1cedc7-8ada-45a5-adef-fbae794dee3e
markercolors = [:green :orange :black :purple :red  :yellow :brown :white]

# ╔═╡ 5fb08873-1ca6-44f2-b68a-438fea6007ed
gr(size=(700,700), xtickfontsize=10, ytickfontsize=10, xguidefontsize=10, yguidefontsize=10, legendfontsize=10, dpi=100, grid=(:y, :gray, :solid, 1, 1.0));


# ╔═╡ 383048a4-12e4-492f-bcbd-a83c6e7dc7aa
Plots.GRBackend()

# ╔═╡ 80542fd1-843e-4d78-9bd3-169c1d6a9672
md"# Notebook"

# ╔═╡ edbcd276-6e0d-41f3-91d8-131b0fc7486b
md"### Read NEMA3 DST"

# ╔═╡ fdff3aeb-834c-4726-8e6c-a35e6b9d46d1
path = datadir("nema3-vac-1m/nema3-vac-1m-0.h5")

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

# ╔═╡ 0c290ab3-9fb0-4b8d-ab2c-c5024f659478
sqrt(pdf.sensor_xyz.x[2]^2+pdf.sensor_xyz.y[2]^2)

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

# ╔═╡ f31a4ae4-11db-4deb-add6-051e8b1db316
begin
	drx = datadir("nema3-vac-1m")
	ffiles = glob("*.h5",drx)
end

# ╔═╡ b76ba213-3a76-4200-bc38-f67ce23f1e5f
DPDF, DPLOR = n3.nema3_analysis(ffiles, 4.0, 1, 10)

# ╔═╡ a9abc4e6-ba71-497f-9370-613a068f4851
dpdf1 = DPDF[1]

# ╔═╡ 77689285-6b75-48f9-ae1e-bbbf710c2bc0
hdf1 = dpdf1["hdf1"]

# ╔═╡ 7d2375e2-9a47-4ca4-8b0d-ca6db6e9a0e8
nrow(hdf1)

# ╔═╡ ed5615f2-e335-4e4d-9b42-306d74f0a9ec
hdf1q = jp.JPetalo.select_by_column_value(hdf1, "q", 4)

# ╔═╡ dd42a228-6d78-4294-a146-f7632c5614e0
nrow(hdf1q)

# ╔═╡ 5233f32d-feea-4ace-ad67-b97681fbc2af
minimum(abs.(hdf1.q))

# ╔═╡ 0513a035-daca-483e-af6e-d03df4c6513d
dx = maximum(abs.(hdf1.x)) - minimum(abs.(hdf1.x))

# ╔═╡ b09c5868-eb8c-478f-a253-eef8aef965d1
dy = maximum(abs.(hdf1.y)) - minimum(abs.(hdf1.y))

# ╔═╡ 361424ff-0d5d-4150-8928-928a184db66f
md"### True info"

# ╔═╡ b943ced5-c805-4dd7-9513-cfe0cb8c47cd
rlors = [dp["rlr"] for dp in DPLOR]

# ╔═╡ 82ea6624-be9b-4619-9efd-2423e117a5f5
function true_hits_from_dpdf(DPDF, tdf="tdf1")
	tth1 = [(dpdf[tdf].event_id[1], dpdf[tdf].x[1], dpdf[tdf].y[1], 
			 dpdf[tdf].z[1], dpdf[tdf].t[1], dpdf[tdf].pre_KE[1]) for dpdf in DPDF]
	event_id = [th[1] for th in tth1]
	x        = [th[2] for th in tth1]
	y        = [th[3] for th in tth1]
	z        = [th[4] for th in tth1]
	t        = [th[5] for th in tth1]
	e        = [th[6] for th in tth1]
	
	return jp.JPetalo.TrueHits(event_id, x, y, z, t, e)
end
	

# ╔═╡ f0c8b379-b0be-494a-a52c-78a3ebec5598
th1 = true_hits_from_dpdf(DPDF,"tdf1")

# ╔═╡ 45a8afe3-1ccb-4b2b-b24e-1a6ad0225446
th2 = true_hits_from_dpdf(DPDF,"tdf2")

# ╔═╡ 2f47d5f9-d1cf-463c-928b-f5b4ca6f6415


# ╔═╡ 5f782529-89e5-43f0-aa74-72614ea184f8
function reco_hits_from_dpdf(DPDF, rdf="hdf1")
	NSIPM = [1]
	DX = [1.0]
	DY = [1.0]
	DZ = [1.0]
	QM = [1]
	for dpf in DPDF
		hdf = dpf[rdf]
		push!(NSIPM, nrow(hdf))
		#append!(Q, hdf.q)
		push!(QM, maximum(hdf.q))
		push!(DX, maximum(abs.(hdf.x)) - minimum(abs.(hdf.x)))
		push!(DY, maximum(abs.(hdf.y)) - minimum(abs.(hdf.y)))
		push!(DZ, maximum(abs.(hdf.z)) - minimum(abs.(hdf.z)))
	end      
	
	return DataFrame(nsipm = NSIPM[2:end], qm = QM[2:end], 
		dx = DX[2:end], dy = DY[2:end], dz = DZ[2:end])
end

# ╔═╡ 3c1e25a5-a4a2-4591-a04c-1e362c7800c5
pxyqth,pxyth,pxzth,pyzth = jp.JPetalo.plot_truehits(th1, th2, 200.0);

# ╔═╡ 3d9e34ea-9c55-45f1-9f41-b159d49f0fa3
begin
plot(pxyqth, layout = (1, 1), aspect_ratio=:equal,size = (1400, 1000), legend=false, fmt = :png)
end

# ╔═╡ 7aaec410-01df-4d32-b919-39f8dca12400
plot(pxyth,pxzth,pyzth, layout = (1, 3), aspect_ratio=:equal,size = (1400, 1000), legend=false, fmt = :png)

# ╔═╡ f5ae929b-b12e-4761-a503-45c5ff750cf2
collect(1:2:10)

# ╔═╡ 4d24d7fa-28f7-4d0a-9889-f0790ae6fe53
tlxy, tlxz, tlyz = jp.JPetalo.plot_lors(th1, th2, 1);

# ╔═╡ 302b25f2-cd53-4918-983d-dc37399507a3
plot(tlxy, tlxz, tlyz, layout = (1, 3), aspect_ratio=:equal,size = (2400, 1000), legend=false, fmt = :png)

# ╔═╡ 178d5a6c-2fb8-49cc-884b-414f86ceaad6
ptsxy, ptsxz, ptsyz = jp.JPetalo.plot_lors_all(th1, th2,101., 10);

# ╔═╡ fa847f97-2a04-49be-bc5f-f6405f274b20
th1

# ╔═╡ 2694fc3f-2883-4af7-ba40-6a5b41ce1147
plot(ptsxy, layout= (1, 1), aspect_ratio=:equal,size = (2400, 1000),legend = false,  fmt = :png)

# ╔═╡ f7b06fae-650c-4ea2-9a12-4a28b775f845
plot(ptsxz, layout= (1, 1), aspect_ratio=:equal, size = (2400, 1000),
	ylims=(-100.0, 100.0), legend = false,  fmt = :png)

# ╔═╡ f1d2767f-c729-4f22-8e83-5cd3161f7d93
plot(ptsyz, layout= (1, 1), size = (1600, 800),
	ylims=(-100.0, 100.0), legend = false,  fmt = :png)

# ╔═╡ 10651565-e21d-4eca-b1b9-80bc87ff691e
md"### Write True Lors to hdf5"

# ╔═╡ 83c574fa-4ed0-4fe8-ad2a-3732db71ccbe
th1

# ╔═╡ de9c8de1-3a59-4945-8595-1615c4a32775
vt12 = collect(zip(th1.t, th2.t, th1.x, th1.y, th1.z, th2.x, th2.y, th2.z))

# ╔═╡ 72acbb63-7be0-45f6-af0d-7d146369c52b
tl = [jp.JPetalo.MlemLor(x...) for x in vt12]

# ╔═╡ 3e9bf527-861e-40a8-bc7d-8dc1639da199
tl

# ╔═╡ 0f5e62a5-d447-4246-82c9-801fe42b4344
jp.JPetalo.write_lors_hdf5(datadir("nema3/truelors2.h5"), tl)

# ╔═╡ f4f30742-121d-49ea-ad14-df9650bf6c6b
md"### Reconstructed data"

# ╔═╡ 9217ad43-35e8-47c1-8227-c1b98a4c9f18
rhits = reco_hits_from_dpdf(DPDF, "hdf1")

# ╔═╡ f8751fb1-d2ac-4801-a935-3710fd7cd94b
begin
pnspim = histogram(rhits.nsipm, bin=20)
ylabel!("frequency")
xlabel!("nof sipm")
end

# ╔═╡ 604beb49-18c0-4a60-8d93-f234be5ed2cd
begin
pqm = histogram(rhits.qm[rhits.qm .< 1000], bin=20)
ylabel!("frequency")
xlabel!("qmax (pes)")
end

# ╔═╡ 7fe45914-9657-4d67-97c0-56b28ac3ff40


# ╔═╡ fc67446f-9d20-486e-9c6f-ea0333a0f150
pdir = datadir("../plots")

# ╔═╡ ffda7937-cdf9-4feb-bd44-72a6d588c86f
savefig(ph1, string(pdir,"/","nsipm.png")) 

# ╔═╡ 87e8e480-51c1-4a8a-b4be-3be87557cffd
begin
phQ = histogram(Q[Q .<100], bins=50, lw=2)
xlabel!("Q (pes)")
end

# ╔═╡ 1a46dd91-8d16-4e20-9050-e5b1ac1e6db4
savefig(phQ, string(pdir,"/","Q.png")) 

# ╔═╡ fb40cce4-89b0-47e3-880e-a62cfa35864b
begin
phQm = histogram(Qmx[Qmx .<1000], bins=30, lw=1)
xlabel!("Qmax (pes)")
end


# ╔═╡ 3aea5d2c-5aa8-4990-85b3-e7f00ad83d03
savefig(phQm, string(pdir,"/","Qmx.png")) 

# ╔═╡ 00137813-b018-47d8-8fc7-d572709d3512
d = Normal(390., 1.2)

# ╔═╡ ef1d4754-48f6-4f11-810e-2be969fa1180
x = rand(d)

# ╔═╡ 1737bcee-e4b6-4582-abfa-536aa9e19cd6
md"## Select data from MC True"

# ╔═╡ e13cd1d9-25b5-4024-a4ff-dc0bacbfd8dc
dr = datadir("nema3-vac-1m")

# ╔═╡ 42679e73-e88f-459b-89ec-401bbf32bb3c
files = glob("*.h5",dr)

# ╔═╡ af7964ec-c92e-46a4-b571-ee75146e1899
typeof(files)

# ╔═╡ 9f44c65e-a304-412e-b22c-80364678c8e9
pdf

# ╔═╡ c589e892-a26d-4e02-9e7a-f4a2772a773d
length(files)

# ╔═╡ 801ded29-0d7f-4407-94e8-7fffb95613b8
th1

# ╔═╡ 8a6c603d-a33b-440b-b5bd-9b7bef37b2b5
md" - First step: select an event"

# ╔═╡ f3b71ecf-fa44-480d-bb7d-612edc6ff15e
#bind evt0 NumberField(0:50000; default=1)

# ╔═╡ 9a42faa0-cbcd-4398-ba7c-2aaf5628e05e
evt0 = 2

# ╔═╡ 7418aa1c-cc4d-47cc-a132-26f01c0aa761
md"- Selected event = $evt0 (use window to change selection)"

# ╔═╡ 4c6aa884-e7f2-4923-84f9-e78869670e1e
evt = jp.JPetalo.select_event(pdf.total_charge, evt0)

# ╔═╡ 9e48d056-0f4e-4e5b-9a8d-3d0e365a798a
sids = evt[!,:sensor_id]

# ╔═╡ 9cc5a574-6cf1-4d1a-ab3a-bd980527940b
pdf.sensor_xyz

# ╔═╡ e2dbad6a-212e-4a32-9b53-1c6aebb13152
pos = jp.JPetalo.sipm_pos.((pdf.sensor_xyz,),sids)

# ╔═╡ 23e66881-21b9-4d8e-92e4-f4540f148e80
length(pos)

# ╔═╡ 970004a5-7ddf-4f9f-97a0-d4ecf105632e
pdf.waveform

# ╔═╡ acabf30d-7cb1-4628-9b9c-4862a2ca86df
function sipm_time(wfm::DataFrame, index::Integer)
	return Array(jp.JPetalo.select_by_index(wfm, "sensor_id", index)[!,3])
end

# ╔═╡ 083e3b1b-7bd1-4654-89d9-9d38109d21b7
sipm_time(pdf.waveform, 1)

# ╔═╡ 6756c8c5-c9da-4ebc-a43e-e92cbf81c9c0
#sipm_time.((pdf.waveform,), sids)

# ╔═╡ da4a8821-401a-4d52-bbe9-01c83db41a16


# ╔═╡ 6feca4b1-9af3-481c-9559-74b604513b07
hitdf = jp.JPetalo.sipm_xyzq(evt, pdf.sensor_xyz)

# ╔═╡ ce216eb0-4352-454d-8099-6dbc0d3b3a55
rr = sqrt.(hitdf.x.^2 + hitdf.y.^2)

# ╔═╡ bb6c863d-1a15-46ea-b835-8a837a16177b
plot(rr, ylims=(395., 396.))

# ╔═╡ df6aaf5c-5cd8-4d2f-a4a1-f367c4164c5b
length(hitdf.q)

# ╔═╡ aa9584a3-9d41-493a-b4a6-0862e3354f80
evt.charge

# ╔═╡ c471cec7-a1cb-4ffe-ac01-8eaf11e2ed1e
md"`plot_xyzq(hdf)` takes a hit DataFrame and returns the scatter plot of the SiPM hits and the tan(y/x) which allows to discriminate one gamma from the other "

# ╔═╡ 7019b621-6ef5-46e8-abb2-f748f005201d
pxyq, pxy, pxz, pyz = jp.JPetalo.plot_xyzq(hitdf, 100.);

# ╔═╡ 5a75154c-bde7-4a2f-af17-63b524ad5958
plot(pxyq,pxy,pxz,pyz, layout = (2, 2), aspect_ratio=:equal,size = (1400, 1000), legend = false,  fmt = :png)


# ╔═╡ 00c704f6-c518-4a05-b49b-b78972231ffe
md"- Histogramming the energy of the SiPMs one can see a large peak at 1 pes. This is the background (diffused light), and suggests a cut at 1 pes"

# ╔═╡ 90777a77-4523-4496-a21f-42338f7e3079
histogram(hitdf.q[hitdf.q.<100], bins=30)

# ╔═╡ 0d8a7dd5-6529-4d0d-a61b-31893cf92262
@bind ecut NumberField(1.0:100.0; default=5.0)

# ╔═╡ b95a86e9-3503-47b3-973a-ee168616b695
typeof(ecut)

# ╔═╡ 5a7d478a-9a79-415d-b99f-d0c548cebeb7
md" - ecut =$ecut pes (can be changed using window)"

# ╔═╡ 4be9f445-e372-4999-a4cb-67565559b6b5
evtQ1 = evt[evt.charge.>ecut,:]

# ╔═╡ ff4ea7b0-dca8-4e8e-b38b-c3963b549695
histogram(evtQ1.charge[evtQ1.charge.<100], bins=30)

# ╔═╡ c82990b5-8e9e-47ec-a80e-4b4067fe97c8


# ╔═╡ 493721a5-7875-48a0-bc46-735021caa292
hitQdf = jp.JPetalo.sipm_xyzq(evtQ1, pdf.sensor_xyz)

# ╔═╡ 5500660e-d169-4fa4-b53f-c32dfdb0e15a
histogram(hitQdf.q, bins=30)

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
b2 = jp.JPetalo.baricenter(hqpdf)

# ╔═╡ 23d0f15e-d54f-4db1-ad35-c4de300db9ba
b1 = jp.JPetalo.baricenter(hqndf)

# ╔═╡ 72b3c6a8-cd90-40de-b1b0-2907588ecc92
md"- Finally we can draw the three proyections of the baricenter, together with the LORs that connect them"

# ╔═╡ 9a07dac6-f038-429c-a51a-a4237532fe82
sxy, syz, sxz = jp.JPetalo.plot_barycenter(bp,bn, pxyQ1,pxzQ1,pyzQ1, 0.1);

# ╔═╡ e7eda113-32ad-47b9-8a74-10f759165e16
plot(sxy,syz,sxz, layout = (1, 3), aspect_ratio=:equal,size = (1400, 1000),legend = false,  fmt = :png)

# ╔═╡ 7941d91e-9807-48b7-b454-f8f192a2695c
md"## Reconstruction of LORS"

# ╔═╡ a4e19c1f-4ccf-47ce-82c5-316004695b62
rlors

# ╔═╡ 231b2f55-eb4e-41f8-8ab4-f528f08ffc28
typeof(rlors)

# ╔═╡ 5b7f3a8c-8dbf-40c4-8a9a-5577d1162293
rlors[1]

# ╔═╡ d1ee23b6-9a32-4da9-b38d-afc2b8131ca6
function plot_lors_all(rlors,  ic=1)
	function plot_xy()
		rl = rlors[1]
	
		pxy  = scatter([rl.x1,rl.x2], [rl.y1,rl.y2],legend=false)

		lsxy = LineSegment([rl.x1,rl.y1],[rl.x2,rl.y2])
		sxy  = plot!(pxy,lsxy)

		for indx in 2:ic:size(rlors)[1]
			rl = rlors[indx]

			pxy = scatter!(pxy, [rl.x1,rl.x2], [rl.y1,rl.y2],legend=false)

			lsxy = LineSegment([rl.x1,rl.y1],[rl.x2,rl.y2])
			sxy = plot!(pxy,lsxy)
		end
		xlabel!("x")
		ylabel!("y")

		return sxy
	end
	function plot_xz()
		rl = rlors[1]

		pxz  = scatter([rl.x1,rl.x2], [rl.z1,rl.z2], legend=false)

		lsxz = LineSegment([rl.x1,rl.z1],[rl.x2,rl.z2])
		sxz  = plot!(pxz,lsxz)

		for indx in 2:ic:size(rlors)[1]
			rl = rlors[indx]

			pxz = scatter!(pxz, [rl.x1,rl.x2], [rl.z1,rl.z2], legend=false)

			lsxz = LineSegment([rl.x1,rl.z1],[rl.x2,rl.z2])
			sxz = plot!(pxz,lsxz)
		end
		xlabel!("x")
		ylabel!("z")

		return sxz
	end

	function plot_yz()
		rl = rlors[1]

		pyz  = scatter([rl.y1,rl.y2], [rl.z1,rl.z2], legend=false)

		lsyz = LineSegment([rl.y1,rl.z1],[rl.y2,rl.z2])
		syz  = plot!(pyz,lsyz)

		for indx in 2:ic:size(rlors)[1]
			rl = rlors[indx]

			pyz = scatter!(pyz, [rl.y1,rl.y2], [rl.z1,rl.z2], legend=false)

			lsyz = LineSegment([rl.y1,rl.z1],[rl.y2,rl.z2])
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


# ╔═╡ b1b81be0-94be-45dc-8135-7501c1813b90
prsxy, prsxz, prsyz = plot_lors_all(rlors,  10)

# ╔═╡ a5ffffe9-5f16-41cf-b5ac-8cfe958d8bd3
plot(prsxy, layout= (1, 1), aspect_ratio=:equal,size = (2400, 1000),legend = false,  fmt = :png)

# ╔═╡ a05044de-4fa5-44a5-ab4f-bacb02a4c5f0
plot(prsxz, layout= (1, 1), aspect_ratio=:equal, size = (2400, 1000),
	ylims=(-100.0, 100.0), legend = false,  fmt = :png)


# ╔═╡ 3d7871cf-6b51-46ad-9a32-77af5a9a803c
plot(prsyz, layout= (1, 1), size = (1600, 800),
	ylims=(-100.0, 100.0), legend = false,  fmt = :png)

# ╔═╡ Cell order:
# ╠═79cfd2fc-9046-11eb-2b13-1b877d57645d
# ╠═d16a3879-acb4-4cce-9b87-faadf4abfece
# ╠═4ee98da9-f3ce-4782-9538-f878f27ed9f7
# ╠═2e5dd476-4e23-427e-8e72-a6af044eb397
# ╠═26a0d201-f36c-4bcf-9219-a2f22c500598
# ╠═ba42ec81-5d26-4759-87e1-2bc004a83ca9
# ╠═ad72864c-0d9f-418e-a646-ba231690b31f
# ╠═08f6e767-a7cb-4fea-af2d-6345019f6c2b
# ╠═3207f446-8643-11eb-37ba-c9aec47fcb8f
# ╠═5115917a-8644-11eb-19fc-0528741ca75d
# ╠═fc8b79a2-8728-11eb-2da7-e3ffa3ceef08
# ╠═ac019a72-dca3-4e7e-9ffe-3cfbede4106e
# ╠═aa9f1b97-21cb-4592-8f47-8549f1305da2
# ╠═0068d015-f292-4bf7-82f9-c6c0115f96e2
# ╠═3ac83772-8efb-4d5f-8639-398d3157e58a
# ╠═68e738e7-88bd-41c2-89e4-594f07d64ddc
# ╠═0f2f4c78-8729-11eb-2bab-27812ce8c47e
# ╠═5ee27d52-86fd-11eb-365e-9f2b1a095575
# ╠═621ec96c-86fd-11eb-1c41-379cc17180dc
# ╠═9b853f27-4288-42a6-8f12-ca004e1773b7
# ╠═26225a47-c7aa-4319-89f2-1abbc5c1116a
# ╠═cf89b973-7b0f-483d-8ce9-ba426f1df2a6
# ╠═d29c1b29-4e5e-44aa-8ca0-334e2a74d6d2
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
# ╠═0c290ab3-9fb0-4b8d-ab2c-c5024f659478
# ╠═fc41144a-ba14-4543-9523-c5450744e126
# ╠═a2f23d32-150b-4bb7-96d0-79f687747e40
# ╠═6fd444d3-03f3-4b56-a02d-cddba85ca96f
# ╠═06208989-8940-4cd5-8f59-c5e7c167b97b
# ╠═5c4b2813-991e-4908-833c-195309dd95f1
# ╠═bfc9e52f-7889-48fc-a6a8-55b8844420e0
# ╠═abfb7e89-99dd-4e2a-8b5f-ec8fefb8fdc3
# ╟─766c0b57-9a16-44d1-ab5f-e7c1d39a5176
# ╠═48a21c9c-cfc3-4d23-b1d9-cb9e47753bb4
# ╠═f31a4ae4-11db-4deb-add6-051e8b1db316
# ╠═b76ba213-3a76-4200-bc38-f67ce23f1e5f
# ╠═a9abc4e6-ba71-497f-9370-613a068f4851
# ╠═77689285-6b75-48f9-ae1e-bbbf710c2bc0
# ╠═7d2375e2-9a47-4ca4-8b0d-ca6db6e9a0e8
# ╠═ed5615f2-e335-4e4d-9b42-306d74f0a9ec
# ╠═dd42a228-6d78-4294-a146-f7632c5614e0
# ╠═5233f32d-feea-4ace-ad67-b97681fbc2af
# ╠═0513a035-daca-483e-af6e-d03df4c6513d
# ╠═b09c5868-eb8c-478f-a253-eef8aef965d1
# ╠═361424ff-0d5d-4150-8928-928a184db66f
# ╠═b943ced5-c805-4dd7-9513-cfe0cb8c47cd
# ╠═82ea6624-be9b-4619-9efd-2423e117a5f5
# ╠═f0c8b379-b0be-494a-a52c-78a3ebec5598
# ╠═45a8afe3-1ccb-4b2b-b24e-1a6ad0225446
# ╠═2f47d5f9-d1cf-463c-928b-f5b4ca6f6415
# ╠═5f782529-89e5-43f0-aa74-72614ea184f8
# ╠═3c1e25a5-a4a2-4591-a04c-1e362c7800c5
# ╠═3d9e34ea-9c55-45f1-9f41-b159d49f0fa3
# ╠═7aaec410-01df-4d32-b919-39f8dca12400
# ╠═f5ae929b-b12e-4761-a503-45c5ff750cf2
# ╠═4d24d7fa-28f7-4d0a-9889-f0790ae6fe53
# ╠═302b25f2-cd53-4918-983d-dc37399507a3
# ╠═178d5a6c-2fb8-49cc-884b-414f86ceaad6
# ╠═fa847f97-2a04-49be-bc5f-f6405f274b20
# ╠═2694fc3f-2883-4af7-ba40-6a5b41ce1147
# ╠═f7b06fae-650c-4ea2-9a12-4a28b775f845
# ╠═f1d2767f-c729-4f22-8e83-5cd3161f7d93
# ╠═10651565-e21d-4eca-b1b9-80bc87ff691e
# ╠═83c574fa-4ed0-4fe8-ad2a-3732db71ccbe
# ╠═de9c8de1-3a59-4945-8595-1615c4a32775
# ╠═72acbb63-7be0-45f6-af0d-7d146369c52b
# ╠═3e9bf527-861e-40a8-bc7d-8dc1639da199
# ╠═0f5e62a5-d447-4246-82c9-801fe42b4344
# ╠═f4f30742-121d-49ea-ad14-df9650bf6c6b
# ╠═9217ad43-35e8-47c1-8227-c1b98a4c9f18
# ╠═f8751fb1-d2ac-4801-a935-3710fd7cd94b
# ╠═604beb49-18c0-4a60-8d93-f234be5ed2cd
# ╠═7fe45914-9657-4d67-97c0-56b28ac3ff40
# ╠═fc67446f-9d20-486e-9c6f-ea0333a0f150
# ╠═ffda7937-cdf9-4feb-bd44-72a6d588c86f
# ╠═87e8e480-51c1-4a8a-b4be-3be87557cffd
# ╠═1a46dd91-8d16-4e20-9050-e5b1ac1e6db4
# ╠═fb40cce4-89b0-47e3-880e-a62cfa35864b
# ╠═3aea5d2c-5aa8-4990-85b3-e7f00ad83d03
# ╠═37dde9a4-3735-441b-b782-478e364468d6
# ╠═c1c9289b-adf3-4987-8b3c-f1eec3b91388
# ╠═d215ef1a-817d-4cc7-8c3a-2ba8d5c411cb
# ╠═00137813-b018-47d8-8fc7-d572709d3512
# ╠═ef1d4754-48f6-4f11-810e-2be969fa1180
# ╠═1737bcee-e4b6-4582-abfa-536aa9e19cd6
# ╠═e13cd1d9-25b5-4024-a4ff-dc0bacbfd8dc
# ╠═42679e73-e88f-459b-89ec-401bbf32bb3c
# ╠═af7964ec-c92e-46a4-b571-ee75146e1899
# ╠═9f44c65e-a304-412e-b22c-80364678c8e9
# ╠═c589e892-a26d-4e02-9e7a-f4a2772a773d
# ╠═801ded29-0d7f-4407-94e8-7fffb95613b8
# ╠═8a6c603d-a33b-440b-b5bd-9b7bef37b2b5
# ╠═f3b71ecf-fa44-480d-bb7d-612edc6ff15e
# ╠═9a42faa0-cbcd-4398-ba7c-2aaf5628e05e
# ╠═7418aa1c-cc4d-47cc-a132-26f01c0aa761
# ╠═4c6aa884-e7f2-4923-84f9-e78869670e1e
# ╠═9e48d056-0f4e-4e5b-9a8d-3d0e365a798a
# ╠═9cc5a574-6cf1-4d1a-ab3a-bd980527940b
# ╠═e2dbad6a-212e-4a32-9b53-1c6aebb13152
# ╠═23e66881-21b9-4d8e-92e4-f4540f148e80
# ╠═970004a5-7ddf-4f9f-97a0-d4ecf105632e
# ╠═acabf30d-7cb1-4628-9b9c-4862a2ca86df
# ╠═083e3b1b-7bd1-4654-89d9-9d38109d21b7
# ╠═6756c8c5-c9da-4ebc-a43e-e92cbf81c9c0
# ╠═da4a8821-401a-4d52-bbe9-01c83db41a16
# ╠═6feca4b1-9af3-481c-9559-74b604513b07
# ╠═ce216eb0-4352-454d-8099-6dbc0d3b3a55
# ╠═bb6c863d-1a15-46ea-b835-8a837a16177b
# ╠═df6aaf5c-5cd8-4d2f-a4a1-f367c4164c5b
# ╠═aa9584a3-9d41-493a-b4a6-0862e3354f80
# ╟─c471cec7-a1cb-4ffe-ac01-8eaf11e2ed1e
# ╠═7019b621-6ef5-46e8-abb2-f748f005201d
# ╠═5a75154c-bde7-4a2f-af17-63b524ad5958
# ╟─00c704f6-c518-4a05-b49b-b78972231ffe
# ╠═90777a77-4523-4496-a21f-42338f7e3079
# ╠═0d8a7dd5-6529-4d0d-a61b-31893cf92262
# ╠═b95a86e9-3503-47b3-973a-ee168616b695
# ╟─5a7d478a-9a79-415d-b99f-d0c548cebeb7
# ╠═4be9f445-e372-4999-a4cb-67565559b6b5
# ╠═ff4ea7b0-dca8-4e8e-b38b-c3963b549695
# ╠═c82990b5-8e9e-47ec-a80e-4b4067fe97c8
# ╠═493721a5-7875-48a0-bc46-735021caa292
# ╠═5500660e-d169-4fa4-b53f-c32dfdb0e15a
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
# ╟─72b3c6a8-cd90-40de-b1b0-2907588ecc92
# ╠═9a07dac6-f038-429c-a51a-a4237532fe82
# ╠═e7eda113-32ad-47b9-8a74-10f759165e16
# ╠═7941d91e-9807-48b7-b454-f8f192a2695c
# ╠═a4e19c1f-4ccf-47ce-82c5-316004695b62
# ╠═231b2f55-eb4e-41f8-8ab4-f528f08ffc28
# ╠═5b7f3a8c-8dbf-40c4-8a9a-5577d1162293
# ╠═d1ee23b6-9a32-4da9-b38d-afc2b8131ca6
# ╠═b1b81be0-94be-45dc-8135-7501c1813b90
# ╠═a5ffffe9-5f16-41cf-b5ac-8cfe958d8bd3
# ╠═a05044de-4fa5-44a5-ab4f-bacb02a4c5f0
# ╠═3d7871cf-6b51-46ad-9a32-77af5a9a803c
