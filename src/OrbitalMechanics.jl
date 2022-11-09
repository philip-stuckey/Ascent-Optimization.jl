module OrbitalMechanics

using StaticArrays
using LinearAlgebra

export Vec, basis

const Vec = SVector{3, Float64}
basis(r⃗) = (normalize(r⃗), [0 1 0; -1 0 0; 0 0 0]*normalize(r⃗) )


export Planet, surface_gravity, surface_speed, circular_orbit_speed
include("planet.jl")


export Ship, 
	setproperty!, 
	mass, 
	mass_flow_rate, 
	thrust, 
	thrust_vector, 
	delta_v, 
	specific_angular_momentum,
	heading
include("ship.jl")

export specific_energy, 
	semi_major_axis, 
	period, 
	eccentricity, 
	eccentricity_vector, 
	apsies, 
	apoapsis, 
	periapsis,
	apsis_velocity, 
	apoapsis_velocity, 
	periapsis_velocity,
	time_since_periapsis,
	time_to_apsis,
	true_anomaly,
	eccentric_anomaly,
	mean_anomaly
include("orbit.jl")


export Simulate!, SimulationParameters, standard_parameters
include("simulation.jl")


export plot_orbit, colTuple, animate_path
include("plotting.jl")

export Maneuver, runManeuver!, burn_time
include("maneuvers.jl")

export Model, runModel
include("MachineLearning/modeling.jl")


end # module OrbitalMechanics
