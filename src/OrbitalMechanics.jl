module OrbitalMechanics

using StaticArrays
using LinearAlgebra

const Vec = SVector{3, Float64}
basis(r⃗) = (normalize(r⃗), [0 1 0; -1 0 0; 0 0 0]*normalize(r⃗) )

include("ship.jl")
include("planet.jl")


function specific_energy(ship::Ship, body::Planet) 
	r = ship.position
	v = ship.velocity
	μ = body.gravitational_parameter
	return (v⋅v)/2 - μ/norm(r)
end



end # module OrbitalMechanics
