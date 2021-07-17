using DataFrames
using LinearAlgebra
using GLM

#linear fit wrapper

function lfit(ndf::DataFrame)
    lr = lm(@formula(y_mean ~ x_mean), ndf)
    c = coef(lr)
    return x -> c[1] + c[2]*x, predict(lr), c
end
#misc math

"""
	rxy(x::Number,y::Number)
	r -> sqrt{x^2 + y^2}
"""
function rxy(x::Number,y::Number)
    return sqrt(x^2 + y^2)
end


"""
	phixy(x::Number,y::Number)
	phi -> atan(y/x)
"""
function phixy(x::Number,y::Number)
    return atan(y,x)
end


"""
	dxyz(x1::Vector{Number}, x2::Vector{Number})

Distance between two points.
"""
function gdxyz(x1, x2)
    return sqrt((x1[1] - x2[1])^2 + (x1[2] - x2[2])^2 + (x1[3] - x2[3])^2)
end


function dxyz(x1::Vector{Float32}, x2::Vector{Float32})
    return gdxyz(x1, x2)
end
function dxyz(x1::Vector{Float64}, x2::Vector{Float64})
    return gdxyz(x1, x2)
end
function dxyz(x1::Vector{Number}, x2::Vector{Number})
    return gdxyz(x1, x2)
end
function dxyz(x1::Vector{Integer}, x2::Vector{Integer})
    return gdxyz(x1, x2)
end


"""
function gwstd(x, q)

Compute the std deviation in x weighted by q:
Sqrt(1/Q Sum_i (x - x_mean) * qi )
"""
function gwstd(x, q)
	xmean = mean(x)
	qs = sum((x.-xmean).^2 .* q)
	Q = sum(q)
	return sqrt(qs/Q)
end
wstd(x::Vector{Float32}, q::Vector{Float32}) = gwstd(x, q)
wstd(x::Vector{Float64}, q::Vector{Float64}) = gwstd(x, q)
wstd(x::Vector{Number}, q::Vector{Number}) = gwstd(x, q)


#q correction
"""
	qcor(df::DataFrame, lf::Function, yref::Number)

Linear correction to q
"""
function qcor(df::DataFrame, lf::Function, yref::Number)
    yold = df.y_mean
    ypre = lf.(df.x_mean)
    return (yold .- ypre) .+ yref
end


function qcor2!(df::DataFrame, lf::Function,
	            xc="r", yc="q1", new="qc", yref=2000.0)
    yold = df[!, yc]
    ypre = lf.(df[!, xc])
    df[!, new] = (yold .- ypre) .+ yref
end


function qcorrection!(n3df::DataFrame, bins::Integer=100)
    prqdf = JPetalo.p1df(n3df.r, n3df.q1, bins)
    lfrq, prq, crq = JPetalo.lfit(prqdf)
    JPetalo.qcor2!(n3df, lfrq, "r", "q1", "qc", 2000.0);
end

#data frames
"""
	select_event(dbdf::DataFrame, index::Int64)

Take the event dataframe and the index of an event and returns a data frame
which selects that particular event

"""
function select_event(dbdf::DataFrame, index::Integer)
	return select_by_index(dbdf, "event_id", index)[:,2:end]
end

"""
	select_by_column_value(df::DataFrame, column::String, value)

Select elements in the DF which have "value" in "column"

"""
function select_by_column_value(df::DataFrame, column::String, value)
	mask = df[!,column].==value
	return df[mask,:]
end

"""
	select_by_column_value_lt(df::DataFrame, column::String, value)

Select elements in the DF which are less than "value" in "column"

"""
function select_by_column_value_lt(df::DataFrame, column::String, value)
	mask = df[!,column].<value
	return df[mask,:]
end


"""
	select_by_column_value_gt(df::DataFrame, column::String, value)

Select elements in the DF which are larger than "value" in "column"

"""
function select_by_column_value_gt(df::DataFrame, column::String, value)
	mask = df[!,column].>value
	return df[mask,:]
end

"""
	select_by_column_value_interval(df::DataFrame, column::String, valuef, valuel)

Select elements in the DF which are in interval (valuef, valuel)

"""
function select_by_column_value_interval(df::DataFrame, column::String, valuef, valuel)
	df1 = select_by_column_value_gt(df, column, valuef)
    return select_by_column_value_lt(df1, column, valuel)
end


"""
	select_by_index(df::DataFrame, column::String, value::Integer)

Select elements in the DF which have "value" (Integer) in "column"

"""
function select_by_index(df::DataFrame, column::String, value::Integer)
	return select_by_column_value(df, column, value)
end


"""
	get_truehits(GP)

Return TrueHits for each one of the gammas in the true lor.

"""
function get_truehits(GP)
	function xyzte(GP)
		eid = [gp.event_id[1] for gp in GP]
		x = [gp.x[1] for gp in GP]
		y = [gp.y[1] for gp in GP]
		z = [gp.z[1] for gp in GP]
		t = [gp.t[1] for gp in GP]
		e = [gp.pre_KE[1] for gp in GP]
		return eid,x,y,z,t,e
	end

	GP1 = []
	GP2 = []
	for gp in GP
		df1 = select_by_column_value(gp, "track_id", 1)
		df2 = select_by_column_value(gp, "track_id", 2)
		push!(GP1,df1)
		push!(GP2,df2)
	end

	return TrueHits(xyzte(GP1)...), TrueHits(xyzte(GP2)...)
end


"""
	sipm_pos(sxyz::DataFrame, index::Integer)

Return the position of the SiPMs in the sipm_xyz database for a given index

"""
function sipm_pos(sxyz::DataFrame, index::Integer)
	return Array(select_by_index(sxyz, "sensor_id", index)[1,2:end])
end


"""
	sipm_xyzq(evt::DataFrame, sxyz::DataFrame)

Return the hits for an event

"""
# function sipm_xyzq(evt::DataFrame, sxyz::DataFrame)
# 	sids = evt[!,:sensor_id]
# 	pos = sipm_pos.((sxyz,),sids)
# 	x = [p[1] for p in pos]
# 	y = [p[2] for p in pos]
# 	z = [p[3] for p in pos]
# 	q = evt[!,:charge]
# 	return DataFrame(x=x,y=y,z=z,q=q)
# end

function sipm_xyzq(qdf::DataFrame, sxyz::DataFrame)
    sids = qdf.sensor_id
    pos = sipm_pos.((sxyz,),sids)
    x = [p[1] for p in pos]
    y = [p[2] for p in pos]
    z = [p[3] for p in pos]
    return DataFrame(x=x,y=y,z=z,q=qdf.Q)
end


"""
	find_xyz_sipm_qmax(hitdf::DataFrame)
Return the coordinates of the SiPM with maximum charge

"""
function find_xyz_sipm_qmax(hitdf::DataFrame)
	qxmax, xmax = find_max_xy(hitdf,"x", "q")
	qymax, ymax = find_max_xy(hitdf,"y", "q")
	qzmax, zmax = find_max_xy(hitdf,"z", "q")
	return Hit(xmax,ymax,zmax,qxmax)
end


""""
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
function find_max_xy(df::DataFrame, xc::String, yc::String)
	ymax, imax = findmax(df[!, yc])
	x_ymax = df[imax, xc]
	return ymax, x_ymax
end

"""
	xyz_dot(hitdf::DataFrame, simax::Hit)

Return the dot product between each SiPM in the event and the SiPM of max charge
"""
function xyz_dot(hitdf::DataFrame, simax::Hit)
	xyzmax = [simax.x, 	simax.y, simax.z]
	xyzm_dot = dot(xyzmax, xyzmax)
	return [dot(Array(hitdf[i,1:3]), xyzmax) /xyzm_dot  for i in 1:nrow(hitdf)]
end


#radius(x::Number, y::Number) = sqrt(x^2 + y^2)
#radius(x::Float64, y::Float64) = sqrt(x^2 + y^2)
