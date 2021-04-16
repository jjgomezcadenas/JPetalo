### A Pluto.jl notebook ###
# v0.14.2

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
import Pkg; Pkg.add.(["Markdown", "Pandas","LazySets"])

# ╔═╡ 3207f446-8643-11eb-37ba-c9aec47fcb8f
begin
	using Markdown
	using InteractiveUtils
	using PlutoUI
	using Test
	using Plots
	using CSV
	using DataFrames
	using Statistics
	using LazySets
	import Pandas
end

# ╔═╡ 5115917a-8644-11eb-19fc-0528741ca75d
begin
	using Unitful
	using UnitfulEquivalences
	using PhysicalConstants.CODATA2018
	using StatsPlots
	using QuadGK
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

# ╔═╡ 5f5d4e88-5819-48a4-b021-4d0a05f7db0f
md"## Types"

# ╔═╡ 9b47dc97-56ef-48ae-8a3b-b5a9a5d5bf6d
"""
Represents a (high level) SiPM hit
"""
struct Hit
	x::Float64
	y::Float64
	z::Float64
	q::Float64
end

# ╔═╡ 56b34074-ac36-476b-b3d9-9057fad68693
md"# Functions"

# ╔═╡ d5299e35-dc5d-4241-91ad-e4f1f7431d06
"""
Selects the elements of a dataframe in terms of value.
"""
function select_by_column_value(df::DataFrame, column::String, value)
	mask = df[!,column].==value
	return df[mask,:]
end

# ╔═╡ e1d24add-2889-48bc-9e05-daf7c9b850aa
"""
Selects the elements of a dataframe in terms of an index
"""
select_by_index(df::DataFrame, column::String, value::Int64) = select_by_column_value(df::DataFrame, column::String, value)

# ╔═╡ 818d753c-cdc9-4489-b7f5-4f32e7f89c4a
"""
Takes the sipm database and the index of a sipm and returns its position
"""
function sipm_pos(dbdf::DataFrame, index::Int64)
	return Array(select_by_index(dbdf, "SensorID", index)[1,2:end])
end

# ╔═╡ b9cc585f-46bd-4f45-9556-e637da68353d
"""
Takes the event dataframe and the index of an event and returns a data frame which selects that particular event

"""
function select_event(dbdf::DataFrame, index::Int64)
	return select_by_index(dbdf, "event_id", index)[:,2:end]
end

# ╔═╡ ed00d3e5-5de3-4446-9d28-3546df3440cf
"""
Takes an event dataframe and a dataframe with the sipm positions and returns a hit-dataframe, with (x,y,z,q) for each sipm

"""
function sipm_xyzq(evt::DataFrame, pdb::DataFrame)
	sids = evt[!,:sensor_id]
	pos = sipm_pos.((pdb,),sids)
	x = [p[1] for p in pos]
	y = [p[2] for p in pos]
	z = [p[3] for p in pos]
	q = evt[!,:charge]
	return DataFrame(x=x,y=y,z=z,q=q)
end

# ╔═╡ 6978a94d-9337-48b8-88f3-7671113747a6
"""
Takes a hit dataframe and splits it into two dataframes depending on whether phi is positve or negative
"""
function fphisel(hdf::DataFrame)
	phi = atan.(hdf.y,hdf.x)
	return hdf[(phi.>0), :], hdf[(phi.<0), :]
end

# ╔═╡ 01e62f30-3713-4e7c-a655-c680fecd7685
"""
Given vectors x and y the function returns the man and std of the positive and negative phi distributions. The types of x and y are nos specified to allow DataFrame rows as well as Arrays. 
"""
function fphi(x,y)
	phi = atan.(y,x)
	phip = phi[phi.>0]
	phim = phi[phi.<0]
	return mean(phip), std(phip), mean(phim), std(phim)
end
	

# ╔═╡ 4274c0db-be70-48e5-b87e-fed080157766
"""
Given a Hit DataFrame (a dataframe with columns x,y,z,q) returns its barycenter
"""
function baricenter(hdf::DataFrame)
	function xq(hdf::DataFrame, pos::String)
		return sum(hdf[!,pos] .* hdf.q) / qt
	end
	qt = sum(hdf.q)
	return Hit(xq(hdf, "x"), xq(hdf, "y"), xq(hdf, "z"), qt)
end

# ╔═╡ d0406ed0-0b88-4568-8dd7-d21d46a5e3b5
md"## Plot functions"

# ╔═╡ 2bd2fa97-51b5-403a-9f76-5282313c7087
"""
Takes a hit dataframe and produces 4 scatter plots showing the clustering of the SiPMs and the phi distribution showing how the clusters distribute at opposite phi values 
"""
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
	pphi = plot(atan.(hdf.y,hdf.x), shape  = :circle, legend=false)
	xlabel!("SiPM number")
	ylabel!("tan(y/x)")
	hphi = histogram(atan.(hdf.y,hdf.x), bins=20)
	xlabel!("tan(y/x)")
	return pxyq,pxy,pxz,pyz, pphi, hphi
end

# ╔═╡ b90b63e2-49da-4003-8fd7-68c81c1df43e
"""
Plots the barycenter of the SiPm clusters together with the LOR that connects them
"""
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

# ╔═╡ 80542fd1-843e-4d78-9bd3-169c1d6a9672
md"# Notebook"

# ╔═╡ edbcd276-6e0d-41f3-91d8-131b0fc7486b
md"### Read Petalo DB"

# ╔═╡ 832c369b-13dd-4003-9381-36c301b5a8de
begin
	dbpath  = datadir("db/petalodb.csv")
	fdf     = CSV.File(dbpath)
	pdb = DataFrame(fdf)[!,2:5]
end

# ╔═╡ d76eada3-4e54-4483-b9f6-321e5740854c
md" - `select_by_index()` on the column SensorId selects a particular sipm position"

# ╔═╡ 689a82e9-0163-4909-b9a2-2a42306f8d6f
select_by_index(pdb, "SensorID", 1000)

# ╔═╡ 606f9c51-b806-494d-93a2-8b6dbbeaad86
md"- We can test that the positions of the sipms corresponding to a given index are correctly obtained"

# ╔═╡ d8023a01-3d19-4741-be7f-41cf13501c55
@test sipm_pos(pdb, 1000) ≈ Array(select_by_index(pdb, "SensorID", 1000)[1,2:end])

# ╔═╡ 8db332bb-9a40-4bdf-b018-64f6a35f6bdf
md"- The scatter plot shows that the positions of the SiPMs form a ring" 

# ╔═╡ b8aacd81-9d8f-4105-a214-220be62b3e43
scatter(pdb.X, pdb.Y, title = "XY positions SiPMs")

# ╔═╡ c0cb1786-d693-42da-b8ca-a2ca4e32cc69
md"- We can alswo show they form a wall in z"

# ╔═╡ 42e1caa1-4aba-44e7-ab32-4538662bcf29
scatter(pdb.X, pdb.Z, title = "XZ positions SiPMs")

# ╔═╡ 4178e1b4-b073-40eb-bda6-c11ecf2b8073
md" ### Read data file"

# ╔═╡ 74e6fc15-9e9f-460f-89ad-f36956b91ae9
begin
	fbp = datadir("fbpet/full_body_phantom_paper.19.csv")
	fbpdf     = CSV.File(fbp)
	snsr = DataFrame(fbpdf)[!,2:end]
end

# ╔═╡ 766c0b57-9a16-44d1-ab5f-e7c1d39a5176
md"- Each data file is indexed by event\_id and sensor\_id. The field sensor_id is used to extract the positions of the sipms in the DB and also provides the charge for a given time bin. This analysis is not considering time (in fact the data is binned with very large data bins, which in practice are irrelevant). Time analysis is illustrated in another notebook"

# ╔═╡ 48a21c9c-cfc3-4d23-b1d9-cb9e47753bb4
md"### Analysis"

# ╔═╡ 8a6c603d-a33b-440b-b5bd-9b7bef37b2b5
md" - First: select an event"

# ╔═╡ f3b71ecf-fa44-480d-bb7d-612edc6ff15e
@bind evt0 NumberField(570000:570100; default=570000)

# ╔═╡ 7418aa1c-cc4d-47cc-a132-26f01c0aa761
md"- Selected event = $evt0 (use window to change selection)"

# ╔═╡ 4c6aa884-e7f2-4923-84f9-e78869670e1e
evt = select_event(snsr, evt0)

# ╔═╡ 9e48d056-0f4e-4e5b-9a8d-3d0e365a798a
sids = evt[!,:sensor_id]

# ╔═╡ 267c0bef-4701-422b-a063-5491780ef016
md" - The number of SiPMs with charge per event is pretty high. For example:
- for event $evt0 
- the number of SiPMs in event =$(length(sids))"

# ╔═╡ a854359d-8442-4d4d-9931-960165f0dd74
md"- Funcion `sipm_pos(sipmdb, index)` takes the name of the SiPM data base (sipmdb) and the index of a sipm and returns the position of the sipm"

# ╔═╡ 4c744ad7-189e-4ed8-96f3-f3e2ba6fa3fa
sipm_pos(pdb, 77194)

# ╔═╡ e2d1b8cf-f2ed-487c-9d30-c957cc9bb0a0
md"- When applied to a vector of indexes the function resturns a vector in which each element is a position"

# ╔═╡ 5570f5de-8a6e-46c2-95ef-0d4b1b172bb0
pos = sipm_pos.((pdb,),sids)

# ╔═╡ ac738ea2-d2fd-4f0b-9b97-dc1745cb8e22
md"Function `sipm_xyzq(evt, pdb)` takes a DataFrame represeting the event and the database (also a DataFrame) to return a hit-DataFrame --a hit is a structure (x,y,z,q)--"

# ╔═╡ 6feca4b1-9af3-481c-9559-74b604513b07
hitdf = sipm_xyzq(evt, pdb)

# ╔═╡ db157d7e-70e5-4d20-8fab-92454b5c2e09
@test hitdf.q ≈ evt.charge

# ╔═╡ 6b778d2c-1b7a-40ca-be5b-c67c1e066e77
begin
	ah = [Array(hitdf[i,1:3]) for i in 1:nrow(hitdf)] 
	@test ah ≈ pos
end

# ╔═╡ c471cec7-a1cb-4ffe-ac01-8eaf11e2ed1e
md"`plot_xyzq(hdf)` takes a hit DataFrame and returns the scatter plot of the SiPM hits and the tan(y/x) which allows to discriminate one gamma from the other "

# ╔═╡ 7019b621-6ef5-46e8-abb2-f748f005201d
pxyq,pxy,pxz,pyz, pphi, hphi = plot_xyzq(hitdf);

# ╔═╡ 5a75154c-bde7-4a2f-af17-63b524ad5958
plot(pxyq,pxy,pxz,pyz,pphi, hphi, layout = (3, 2), legend = false,  fmt = :png)


# ╔═╡ d7f41a07-7c9e-43ae-8e35-48bc62bb8c08
md"- Notice the large dispersion in hits, in particular in the phi distributions. This is reflected by the std of the positive/negative phis"

# ╔═╡ 8ce50ad3-318b-4c11-8ec5-753d33cc8451
fphi(hitdf.x, hitdf.y)

# ╔═╡ 00c704f6-c518-4a05-b49b-b78972231ffe
md"- Histogramming the energy of the SiPMs one can see a large peak at 1 pes. This is the background (diffused light), and suggests a cut at 1 pes"

# ╔═╡ 90777a77-4523-4496-a21f-42338f7e3079
histogram(hitdf.q[hitdf.q.<10], bins=30)

# ╔═╡ 0d8a7dd5-6529-4d0d-a61b-31893cf92262
@bind ecut NumberField(1:3; default=1)

# ╔═╡ 5a7d478a-9a79-415d-b99f-d0c548cebeb7
md" - ecut =$ecut pes (can be changed using window)"

# ╔═╡ 4be9f445-e372-4999-a4cb-67565559b6b5
evtQ1 = evt[evt.charge.>ecut,:]

# ╔═╡ 785cd6a8-fd02-4624-b0e8-b487d3186800
hitQdf = sipm_xyzq(evtQ1,pdb)

# ╔═╡ b7625433-4509-422b-9cdd-83b9984f43da
md"- Repeating the plots ater cut at $ecut pes shows a much cleaner distribution"

# ╔═╡ 9f59a9cc-b963-4bd2-b479-e3cdbaa53d8f
pxyqQ1,pxyQ1,pxzQ1,pyzQ1, pphiQ1, hphiQ1 = plot_xyzq(hitQdf);

# ╔═╡ d7d8f4af-bb0e-42c3-aa5d-c657bb17d2af

plot(pxyqQ1,pxyQ1,pxzQ1,pyzQ1,pphiQ1, hphiQ1, layout = (3, 2), legend = false, fmt = :png)

# ╔═╡ cd3f4960-2cbd-4791-8fb6-86d3671d9937
md"- This can also be seen by the much lower value of the phi distributions std"

# ╔═╡ d89bc875-7d74-4ad9-ad7f-1b8f100ee3e3
fphi(hitQdf.x, hitQdf.y)

# ╔═╡ 35132119-c793-4bfe-b228-7a017ce7789d
md"- We can now select the hits with positive and negative phi, which define the clusters of the gammas"

# ╔═╡ 2af988f3-4754-4de8-833b-a5bc57f0381d
hqpdf, hqndf = fphisel(hitQdf);

# ╔═╡ da5306f6-9bea-4cc3-9bf9-767b0908fb69
md"- And compute the baricenters, which define the gamma position"

# ╔═╡ 7f277fff-efaa-4502-8a58-45c4d4a514f5
bp = baricenter(hqpdf)

# ╔═╡ dbd61080-79d6-4e64-8a87-a1ed243ff503
bn = baricenter(hqndf)

# ╔═╡ 72b3c6a8-cd90-40de-b1b0-2907588ecc92
md"- Finally we can draw the three proyections of the baricenter, together with the LORs that connect them"

# ╔═╡ 9a07dac6-f038-429c-a51a-a4237532fe82
sxy, syz, sxz = plot_barycenter(bp,bn, pxy, pxz, pyz);

# ╔═╡ e7eda113-32ad-47b9-8a74-10f759165e16
plot(sxy,syz,sxz, layout = (1, 3), legend = false,  fmt = :png)

# ╔═╡ Cell order:
# ╠═79cfd2fc-9046-11eb-2b13-1b877d57645d
# ╠═d16a3879-acb4-4cce-9b87-faadf4abfece
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
# ╟─5f5d4e88-5819-48a4-b021-4d0a05f7db0f
# ╠═9b47dc97-56ef-48ae-8a3b-b5a9a5d5bf6d
# ╠═56b34074-ac36-476b-b3d9-9057fad68693
# ╠═d5299e35-dc5d-4241-91ad-e4f1f7431d06
# ╠═e1d24add-2889-48bc-9e05-daf7c9b850aa
# ╠═818d753c-cdc9-4489-b7f5-4f32e7f89c4a
# ╠═b9cc585f-46bd-4f45-9556-e637da68353d
# ╠═ed00d3e5-5de3-4446-9d28-3546df3440cf
# ╠═6978a94d-9337-48b8-88f3-7671113747a6
# ╠═01e62f30-3713-4e7c-a655-c680fecd7685
# ╠═4274c0db-be70-48e5-b87e-fed080157766
# ╟─d0406ed0-0b88-4568-8dd7-d21d46a5e3b5
# ╠═2bd2fa97-51b5-403a-9f76-5282313c7087
# ╠═b90b63e2-49da-4003-8fd7-68c81c1df43e
# ╟─80542fd1-843e-4d78-9bd3-169c1d6a9672
# ╟─edbcd276-6e0d-41f3-91d8-131b0fc7486b
# ╠═832c369b-13dd-4003-9381-36c301b5a8de
# ╟─d76eada3-4e54-4483-b9f6-321e5740854c
# ╠═689a82e9-0163-4909-b9a2-2a42306f8d6f
# ╟─606f9c51-b806-494d-93a2-8b6dbbeaad86
# ╠═d8023a01-3d19-4741-be7f-41cf13501c55
# ╟─8db332bb-9a40-4bdf-b018-64f6a35f6bdf
# ╠═b8aacd81-9d8f-4105-a214-220be62b3e43
# ╟─c0cb1786-d693-42da-b8ca-a2ca4e32cc69
# ╠═42e1caa1-4aba-44e7-ab32-4538662bcf29
# ╟─4178e1b4-b073-40eb-bda6-c11ecf2b8073
# ╠═74e6fc15-9e9f-460f-89ad-f36956b91ae9
# ╟─766c0b57-9a16-44d1-ab5f-e7c1d39a5176
# ╟─48a21c9c-cfc3-4d23-b1d9-cb9e47753bb4
# ╟─8a6c603d-a33b-440b-b5bd-9b7bef37b2b5
# ╟─f3b71ecf-fa44-480d-bb7d-612edc6ff15e
# ╟─7418aa1c-cc4d-47cc-a132-26f01c0aa761
# ╠═4c6aa884-e7f2-4923-84f9-e78869670e1e
# ╠═9e48d056-0f4e-4e5b-9a8d-3d0e365a798a
# ╟─267c0bef-4701-422b-a063-5491780ef016
# ╟─a854359d-8442-4d4d-9931-960165f0dd74
# ╠═4c744ad7-189e-4ed8-96f3-f3e2ba6fa3fa
# ╟─e2d1b8cf-f2ed-487c-9d30-c957cc9bb0a0
# ╠═5570f5de-8a6e-46c2-95ef-0d4b1b172bb0
# ╟─ac738ea2-d2fd-4f0b-9b97-dc1745cb8e22
# ╠═6feca4b1-9af3-481c-9559-74b604513b07
# ╠═db157d7e-70e5-4d20-8fab-92454b5c2e09
# ╠═6b778d2c-1b7a-40ca-be5b-c67c1e066e77
# ╟─c471cec7-a1cb-4ffe-ac01-8eaf11e2ed1e
# ╠═7019b621-6ef5-46e8-abb2-f748f005201d
# ╠═5a75154c-bde7-4a2f-af17-63b524ad5958
# ╟─d7f41a07-7c9e-43ae-8e35-48bc62bb8c08
# ╠═8ce50ad3-318b-4c11-8ec5-753d33cc8451
# ╟─00c704f6-c518-4a05-b49b-b78972231ffe
# ╠═90777a77-4523-4496-a21f-42338f7e3079
# ╟─0d8a7dd5-6529-4d0d-a61b-31893cf92262
# ╟─5a7d478a-9a79-415d-b99f-d0c548cebeb7
# ╠═4be9f445-e372-4999-a4cb-67565559b6b5
# ╠═785cd6a8-fd02-4624-b0e8-b487d3186800
# ╠═b7625433-4509-422b-9cdd-83b9984f43da
# ╠═9f59a9cc-b963-4bd2-b479-e3cdbaa53d8f
# ╠═d7d8f4af-bb0e-42c3-aa5d-c657bb17d2af
# ╟─cd3f4960-2cbd-4791-8fb6-86d3671d9937
# ╠═d89bc875-7d74-4ad9-ad7f-1b8f100ee3e3
# ╟─35132119-c793-4bfe-b228-7a017ce7789d
# ╠═2af988f3-4754-4de8-833b-a5bc57f0381d
# ╟─da5306f6-9bea-4cc3-9bf9-767b0908fb69
# ╠═7f277fff-efaa-4502-8a58-45c4d4a514f5
# ╠═dbd61080-79d6-4e64-8a87-a1ed243ff503
# ╟─72b3c6a8-cd90-40de-b1b0-2907588ecc92
# ╠═9a07dac6-f038-429c-a51a-a4237532fe82
# ╠═e7eda113-32ad-47b9-8a74-10f759165e16
