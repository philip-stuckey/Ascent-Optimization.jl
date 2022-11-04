module OrbitalMechanics

using StaticArrays
using LinearAlgebra

export Vec, basis

const Vec = SVector{3, Float64}
basis(r⃗) = (normalize(r⃗), [0 1 0; -1 0 0; 0 0 0]*normalize(r⃗) )


export Ship, 
	setproperty!, 
	mass, 
	mass_flow_rate, 
	thrust, 
	thrust_vector, 
	delta_v, 
	specific_angular_momentum
include("ship.jl")


export Planet, surface_gravity, surface_speed, circular_orbit_speed
include("planet.jl")

export specific_energy, 
	semi_major_axis, 
	period, 
	eccentricity, 
	eccentricity_vector, 
	apsies, 
	apoapsis, 
	periapsis
include("orbit.jl")


export Simulate!
include("simulation.jl")


export plot_orbit
include("plotting.jl")

end # module OrbitalMechanics
