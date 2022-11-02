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

function semi_major_axis(ship,body) 
	μ = body.gravitational_parameter
	ε = specific_energy(ship,body)
	return -μ/(2*ε)
end

function eccentricity(ship,body) 
	μ = body.gravitational_parameter
	ε = specific_energy(ship,body)
	l = specific_angular_momentum(ship)
	return √(1 + 2ε*(l⋅l)/μ^2)
end

function eccentricity_vector(ship,body)
	r = ship.position
	r̂ = normalize(ship.position)
	v = ship.velocity
	l = specific_angular_momentum(ship)
	μ = body.gravitational_parameter
	return r̂ - (v×l)/μ
end


±(a,b) = (a+b, a-b)
function apsies(ship; body::Planet)
    μ = body.gravitational_parameter
    ε = specific_energy(ship,body)
    l = specific_angular_momentum(ship)
    return (-μ ± √(μ^2 + 2ε*(l⋅l))) ./ 2ε
end

apoapsis(ship,body) = maximum(apsies(ship;body))  
periapsis(ship,body) = minimum(apsies(ship;body))


end # module OrbitalMechanics
