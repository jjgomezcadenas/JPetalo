#TYPES

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
