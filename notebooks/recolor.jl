### A Pluto.jl notebook ###
# v0.14.4

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

# ╔═╡ 79cfd2fc-9046-11eb-2b13-1b877d57645d
md"# RECO LOR

- Reconstructs the two gammas in the event and draws the Line Of Response (LOR)
"

# ╔═╡ 4ee98da9-f3ce-4782-9538-f878f27ed9f7


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

# ╔═╡ cf89b973-7b0f-483d-8ce9-ba426f1df2a6
#fbi = ingredients(srcdir("fbi.jl"))

# ╔═╡ 0b1cedc7-8ada-45a5-adef-fbae794dee3e
markercolors = [:green :orange :black :purple :red  :yellow :brown :white]

# ╔═╡ 5fb08873-1ca6-44f2-b68a-438fea6007ed
gr(size=(700,700), xtickfontsize=12, ytickfontsize=12, xguidefontsize=14, yguidefontsize=14, legendfontsize=10, dpi=100, grid=(:y, :gray, :solid, 2, 0.4));


# ╔═╡ 80542fd1-843e-4d78-9bd3-169c1d6a9672
md"# Notebook"

# ╔═╡ edbcd276-6e0d-41f3-91d8-131b0fc7486b
md"### Read Petalo DB"

# ╔═╡ bb32afef-b829-4222-bd12-0372cafa2ccb
@bind option Select(["readFromWeb", "readFromFile"])

# ╔═╡ 61845b04-e1b4-4cdb-be9b-9ee9ec53bb69
if option == "readFromWeb"
	res = HTTP.get("http://next.ific.uv.es/icgallery/public/PETALO/petalodb.csv")
	fdf = CSV.File(res.body)
else
	dbpath  = datadir("db/petalodb.csv")
	fdf     = CSV.File(dbpath)
end

# ╔═╡ 832c369b-13dd-4003-9381-36c301b5a8de
pdb = DataFrame(fdf)[!,2:5]

# ╔═╡ 8db332bb-9a40-4bdf-b018-64f6a35f6bdf
md"- The scatter plot shows that the positions of the SiPMs form a ring" 

# ╔═╡ b8aacd81-9d8f-4105-a214-220be62b3e43
scatter(pdb.X, pdb.Y, title = "XY positions SiPMs", leg=false)

# ╔═╡ c0cb1786-d693-42da-b8ca-a2ca4e32cc69
md"- We can alswo show they form a wall in z"

# ╔═╡ 42e1caa1-4aba-44e7-ab32-4538662bcf29
scatter(pdb.X, pdb.Z, title = "XZ positions SiPMs", leg=false)

# ╔═╡ 4178e1b4-b073-40eb-bda6-c11ecf2b8073
md" ### Read data file"

# ╔═╡ c661c09e-2b68-4bf5-9793-5f2cbb31c025
if option == "readFromWeb"
	res2 = HTTP.get("http://next.ific.uv.es/icgallery/public/PETALO/full_body_phantom_paper.19.h5.csv")
	fbpdf = CSV.File(res2.body)
else
	fbp = datadir("fbpet/full_body_phantom_paper.19.h5.csv")
	fbpdf     = CSV.File(fbp)
end

# ╔═╡ 74e6fc15-9e9f-460f-89ad-f36956b91ae9
snsr = DataFrame(fbpdf)[!,2:end]

# ╔═╡ 766c0b57-9a16-44d1-ab5f-e7c1d39a5176
md"- Each data file is indexed by event\_id and sensor\_id. The field sensor_id is used to extract the positions of the sipms in the DB and also provides the charge for a given time bin. This analysis is not considering time (in fact the data is binned with very large data bins, which in practice are irrelevant). Time analysis is illustrated in another notebook"

# ╔═╡ 48a21c9c-cfc3-4d23-b1d9-cb9e47753bb4
md"### Analysis"

# ╔═╡ 8a6c603d-a33b-440b-b5bd-9b7bef37b2b5
md" - First step: select an event"

# ╔═╡ f3b71ecf-fa44-480d-bb7d-612edc6ff15e
@bind evt0 NumberField(570000:570100; default=570000)

# ╔═╡ 7418aa1c-cc4d-47cc-a132-26f01c0aa761
md"- Selected event = $evt0 (use window to change selection)"

# ╔═╡ ac738ea2-d2fd-4f0b-9b97-dc1745cb8e22
md"Function `sipm_xyzq(evt, pdb)` takes a DataFrame represeting the event and the database (also a DataFrame) to return a hit-DataFrame --a hit is a structure (x,y,z,q)--"

# ╔═╡ c471cec7-a1cb-4ffe-ac01-8eaf11e2ed1e
md"`plot_xyzq(hdf)` takes a hit DataFrame and returns the scatter plot of the SiPM hits and the tan(y/x) which allows to discriminate one gamma from the other "

# ╔═╡ d7f41a07-7c9e-43ae-8e35-48bc62bb8c08
md"- Notice the large dispersion in hits, in particular in the phi distributions. This is reflected by the std of the positive/negative phis"

# ╔═╡ 00c704f6-c518-4a05-b49b-b78972231ffe
md"- Histogramming the energy of the SiPMs one can see a large peak at 1 pes. This is the background (diffused light), and suggests a cut at 1 pes"

# ╔═╡ 0d8a7dd5-6529-4d0d-a61b-31893cf92262
@bind ecut NumberField(1:3; default=1)

# ╔═╡ 5a7d478a-9a79-415d-b99f-d0c548cebeb7
md" - ecut =$ecut pes (can be changed using window)"

# ╔═╡ 493721a5-7875-48a0-bc46-735021caa292


# ╔═╡ b7625433-4509-422b-9cdd-83b9984f43da
md"- Repeating the plots ater cut at $ecut pes shows a much cleaner distribution"

# ╔═╡ cd3f4960-2cbd-4791-8fb6-86d3671d9937
md"- This can also be seen by the much lower value of the phi distributions std"

# ╔═╡ 35132119-c793-4bfe-b228-7a017ce7789d
md"- We can now select the hits with positive and negative phi, which define the clusters of the gammas"

# ╔═╡ da5306f6-9bea-4cc3-9bf9-767b0908fb69
md"- And compute the baricenters, which define the gamma position"

# ╔═╡ 72b3c6a8-cd90-40de-b1b0-2907588ecc92
md"- Finally we can draw the three proyections of the baricenter, together with the LORs that connect them"

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

# ╔═╡ 982329e8-61df-4bb9-b5fa-82f5c50edafe
function test_find_max_xy(evt)
	qmax, iqmax = findmax(evt.charge)
	simax = evt.sensor_id[iqmax]
	qm, sim = find_max_xy(evt,"sensor_id", "charge")
	return qm ≈ qmax && sim ≈ simax
end

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

# ╔═╡ 1cdcefbf-99e1-43b3-9587-b9b5bfee40c4
md" - `simax` is the Hit with the coordinates and charge of the SiPM of max charge"

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

# ╔═╡ 4c6aa884-e7f2-4923-84f9-e78869670e1e
evt = select_event(snsr, evt0)

# ╔═╡ 9e48d056-0f4e-4e5b-9a8d-3d0e365a798a
sids = evt[!,:sensor_id]

# ╔═╡ aa9584a3-9d41-493a-b4a6-0862e3354f80
evt.charge

# ╔═╡ 4be9f445-e372-4999-a4cb-67565559b6b5
evtQ1 = evt[evt.charge.>ecut,:]

# ╔═╡ 6e112f16-c32c-491a-9b68-180a57aace8f
find_max_xy(evt,"sensor_id", "charge")

# ╔═╡ c91baad0-7c07-4b40-a99e-3c773b18b065
scatter(evt.sensor_id, evt.charge, marker=:circle, leg=false)

# ╔═╡ 3f26bf3b-a006-4956-9674-09b157c247c9
@test test_find_max_xy(evt)

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

# ╔═╡ 6feca4b1-9af3-481c-9559-74b604513b07
hitdf = sipm_xyzq(evt, pdb)

# ╔═╡ df6aaf5c-5cd8-4d2f-a4a1-f367c4164c5b
length(hitdf.q)

# ╔═╡ db157d7e-70e5-4d20-8fab-92454b5c2e09
@test hitdf.q == evt.charge

# ╔═╡ 90777a77-4523-4496-a21f-42338f7e3079
histogram(hitdf.q[hitdf.q.<10], bins=30)

# ╔═╡ 44af890b-5d5b-4e57-b284-0a7729674466
test_find_xyz_sipm_qmax(hitdf)

# ╔═╡ e694f819-3099-4a90-9b44-e2151112941c
simax = find_xyz_sipm_qmax(hitdf::DataFrame)

# ╔═╡ 68847964-7edf-4dba-9cd0-2e9197fdf5a9
npr = xyz_dot(hitdf, simax)

# ╔═╡ bbe9358e-1d79-4bd2-b003-d49831e21cfc
histogram(npr)

# ╔═╡ 072d42e4-9924-4d2d-8302-c2cd34ec9ee3
@test maximum(npr) < π/2

# ╔═╡ 2c52dc80-04ca-48e3-acfc-3d6f63acb323
@test minimum(npr) > -π/2 

# ╔═╡ 33192e19-2a20-4a43-ba36-2a5d40d009cf
hitp, hitn  = sipmsel(hitdf)

# ╔═╡ 12044636-6967-424f-977f-c11e8df0b2e1
scatter(fphi(hitdf), leg=false)

# ╔═╡ 785cd6a8-fd02-4624-b0e8-b487d3186800
hitQdf = sipm_xyzq(evtQ1,pdb)

# ╔═╡ 2af988f3-4754-4de8-833b-a5bc57f0381d
hqpdf, hqndf = sipmsel(hitQdf);

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

# ╔═╡ 7f277fff-efaa-4502-8a58-45c4d4a514f5
bp = baricenter(hqpdf)

# ╔═╡ dbd61080-79d6-4e64-8a87-a1ed243ff503
bn = baricenter(hqndf)

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

# ╔═╡ 7019b621-6ef5-46e8-abb2-f748f005201d
pxyq,pxy,pxz,pyz, pphi, hphi = plot_xyzq(hitdf);

# ╔═╡ 5a75154c-bde7-4a2f-af17-63b524ad5958
plot(pxyq,pxy,pxz,pyz,pphi, hphi, layout = (3, 2), legend = false,  fmt = :png)


# ╔═╡ 9f59a9cc-b963-4bd2-b479-e3cdbaa53d8f
pxyqQ1,pxyQ1,pxzQ1,pyzQ1, pphiQ1, hphiQ1 = plot_xyzq(hitQdf);

# ╔═╡ d7d8f4af-bb0e-42c3-aa5d-c657bb17d2af

plot(pxyqQ1,pxyQ1,pxzQ1,pyzQ1,pphiQ1, hphiQ1, layout = (3, 2), legend = false, fmt = :png)

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

# ╔═╡ 9a07dac6-f038-429c-a51a-a4237532fe82
sxy, syz, sxz = plot_barycenter(bp,bn, pxy, pxz, pyz);

# ╔═╡ e7eda113-32ad-47b9-8a74-10f759165e16
plot(sxy,syz,sxz, layout = (1, 3), legend = false,  fmt = :png)

# ╔═╡ Cell order:
# ╠═79cfd2fc-9046-11eb-2b13-1b877d57645d
# ╠═d16a3879-acb4-4cce-9b87-faadf4abfece
# ╠═4ee98da9-f3ce-4782-9538-f878f27ed9f7
# ╠═3207f446-8643-11eb-37ba-c9aec47fcb8f
# ╠═5115917a-8644-11eb-19fc-0528741ca75d
# ╠═fc8b79a2-8728-11eb-2da7-e3ffa3ceef08
# ╠═68e738e7-88bd-41c2-89e4-594f07d64ddc
# ╠═0f2f4c78-8729-11eb-2bab-27812ce8c47e
# ╠═5ee27d52-86fd-11eb-365e-9f2b1a095575
# ╠═621ec96c-86fd-11eb-1c41-379cc17180dc
# ╠═9b853f27-4288-42a6-8f12-ca004e1773b7
# ╠═cf89b973-7b0f-483d-8ce9-ba426f1df2a6
# ╠═0b1cedc7-8ada-45a5-adef-fbae794dee3e
# ╠═5fb08873-1ca6-44f2-b68a-438fea6007ed
# ╟─80542fd1-843e-4d78-9bd3-169c1d6a9672
# ╟─edbcd276-6e0d-41f3-91d8-131b0fc7486b
# ╟─bb32afef-b829-4222-bd12-0372cafa2ccb
# ╠═61845b04-e1b4-4cdb-be9b-9ee9ec53bb69
# ╠═832c369b-13dd-4003-9381-36c301b5a8de
# ╟─8db332bb-9a40-4bdf-b018-64f6a35f6bdf
# ╠═b8aacd81-9d8f-4105-a214-220be62b3e43
# ╟─c0cb1786-d693-42da-b8ca-a2ca4e32cc69
# ╠═42e1caa1-4aba-44e7-ab32-4538662bcf29
# ╟─4178e1b4-b073-40eb-bda6-c11ecf2b8073
# ╠═c661c09e-2b68-4bf5-9793-5f2cbb31c025
# ╠═74e6fc15-9e9f-460f-89ad-f36956b91ae9
# ╟─766c0b57-9a16-44d1-ab5f-e7c1d39a5176
# ╟─48a21c9c-cfc3-4d23-b1d9-cb9e47753bb4
# ╠═8a6c603d-a33b-440b-b5bd-9b7bef37b2b5
# ╠═f3b71ecf-fa44-480d-bb7d-612edc6ff15e
# ╟─7418aa1c-cc4d-47cc-a132-26f01c0aa761
# ╠═4c6aa884-e7f2-4923-84f9-e78869670e1e
# ╠═9e48d056-0f4e-4e5b-9a8d-3d0e365a798a
# ╟─ac738ea2-d2fd-4f0b-9b97-dc1745cb8e22
# ╠═6feca4b1-9af3-481c-9559-74b604513b07
# ╠═df6aaf5c-5cd8-4d2f-a4a1-f367c4164c5b
# ╠═aa9584a3-9d41-493a-b4a6-0862e3354f80
# ╠═db157d7e-70e5-4d20-8fab-92454b5c2e09
# ╟─c471cec7-a1cb-4ffe-ac01-8eaf11e2ed1e
# ╠═7019b621-6ef5-46e8-abb2-f748f005201d
# ╠═5a75154c-bde7-4a2f-af17-63b524ad5958
# ╟─d7f41a07-7c9e-43ae-8e35-48bc62bb8c08
# ╟─00c704f6-c518-4a05-b49b-b78972231ffe
# ╠═90777a77-4523-4496-a21f-42338f7e3079
# ╟─0d8a7dd5-6529-4d0d-a61b-31893cf92262
# ╟─5a7d478a-9a79-415d-b99f-d0c548cebeb7
# ╠═4be9f445-e372-4999-a4cb-67565559b6b5
# ╠═493721a5-7875-48a0-bc46-735021caa292
# ╠═785cd6a8-fd02-4624-b0e8-b487d3186800
# ╟─b7625433-4509-422b-9cdd-83b9984f43da
# ╠═9f59a9cc-b963-4bd2-b479-e3cdbaa53d8f
# ╠═d7d8f4af-bb0e-42c3-aa5d-c657bb17d2af
# ╟─cd3f4960-2cbd-4791-8fb6-86d3671d9937
# ╟─35132119-c793-4bfe-b228-7a017ce7789d
# ╠═2af988f3-4754-4de8-833b-a5bc57f0381d
# ╟─da5306f6-9bea-4cc3-9bf9-767b0908fb69
# ╠═7f277fff-efaa-4502-8a58-45c4d4a514f5
# ╠═dbd61080-79d6-4e64-8a87-a1ed243ff503
# ╟─72b3c6a8-cd90-40de-b1b0-2907588ecc92
# ╠═9a07dac6-f038-429c-a51a-a4237532fe82
# ╠═e7eda113-32ad-47b9-8a74-10f759165e16
# ╟─1709ae5d-b256-4ab6-8c16-5edaa354f867
# ╟─ccf11f8d-5bcb-40e1-8b6d-d88812c41f88
# ╠═35d0f27d-0f97-401a-9a82-9e176cf62fa2
# ╟─56b34074-ac36-476b-b3d9-9057fad68693
# ╟─a124b244-f873-4847-b5de-a18b80975e80
# ╠═b7a9953e-4392-4762-a6f4-979d47426639
# ╟─567ca9fc-ce60-4b30-9246-07d999cd7654
# ╠═12e87e3c-5521-4d94-b7f6-009abf5e96f1
# ╠═6e112f16-c32c-491a-9b68-180a57aace8f
# ╠═c91baad0-7c07-4b40-a99e-3c773b18b065
# ╠═982329e8-61df-4bb9-b5fa-82f5c50edafe
# ╠═3f26bf3b-a006-4956-9674-09b157c247c9
# ╟─300b1614-741d-4720-b2bc-01539423a163
# ╠═3872788b-62a9-4733-8d88-b4af0e59c423
# ╠═d09f2c56-4ac1-47b4-b5d2-534a4c7b5301
# ╠═44af890b-5d5b-4e57-b284-0a7729674466
# ╟─1cdcefbf-99e1-43b3-9587-b9b5bfee40c4
# ╠═e694f819-3099-4a90-9b44-e2151112941c
# ╟─0a2569da-b161-428b-b092-722e01255971
# ╠═735c6816-a7fc-457f-ae8b-95886fcc84bb
# ╟─81685b61-d54e-4c6b-a087-c53a7d0f62e3
# ╠═68847964-7edf-4dba-9cd0-2e9197fdf5a9
# ╠═bbe9358e-1d79-4bd2-b003-d49831e21cfc
# ╠═072d42e4-9924-4d2d-8302-c2cd34ec9ee3
# ╠═2c52dc80-04ca-48e3-acfc-3d6f63acb323
# ╟─8d38fc77-7801-42b6-91d4-d56e3f657f1b
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
# ╟─138db3d0-4920-4405-949f-d1ff3308c696
# ╠═dfcfd7e7-b6f8-423c-a177-be6360a031c1
# ╟─86078e18-ccbf-40cd-9b4c-08f1a5469147
# ╠═dc324771-35b0-4410-a0d5-8582d8975d3a
# ╟─7623540b-7892-45b0-bd28-39c2fedf620d
# ╠═f743a6f0-7f28-4bb7-be72-b25424879312
# ╟─d859296b-b081-4e3b-aa9c-08c8f6fc99cd
# ╟─02434a19-c1cf-4c60-9457-928d64e3419b
# ╠═625b20ac-1c95-4206-85bb-14705a4028f4
# ╟─afac8cb0-83d3-4be8-90a4-329b15d1b614
# ╠═856e70b5-0f27-45c6-afa2-b86cfd19f1bb
