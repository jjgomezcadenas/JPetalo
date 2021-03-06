#TYPES
using PhysicalConstants
using PhysicalConstants.CODATA2018
using Unitful

import Unitful:
    nm, μm, mm, cm, m, km, inch, ft, mi,
    ac,
    mg, g, kg,
    Ra, °F, °C, K,
    rad, °,
    ns, μs, ms, ps, s, minute, hr, d, yr, Hz,
    eV,
    μJ, mJ, J,
	mW, μW, W,
    A, N, mol, mmol, V, L, mL, μL, M

"""
	TrueHits

Represents vectors of true hits, characterized by position and energy

"""
struct TrueHits
	event_id::Vector{Integer}
	x::Vector{Float32}
	y::Vector{Float32}
	z::Vector{Float32}
	t::Vector{Float32}
	e::Vector{Float32}
end

"""
	TrueHit

Represents a true hit, characterized by position and energy

"""
struct TrueHit
	event_id::Integer
	x::Float32
	y::Float32
	z::Float32
	t::Float32
	e::Float32
end

"""
	Hit

Represents a hit, characterized by position and charge

"""
struct Hit
	x::Float32
	y::Float32
	z::Float32
	q::Float32
end

"""
	DetConf

Defines detector configuration
"""
struct DetConf
	pde::Float32
	sigma_tof::Float32
	ecut::Float32
	qmin::Float32
	qmax::Float32
	max_pes::Integer
	ntof::Integer
end
