module AscentOptimization

using StaticArrays
using LinearAlgebra

export Vec, basis

include("units.jl")

const Vec{T} = SVector{3, T}
const rotation = SMatrix{3,3}([0 1 0; -1 0 0; 0 0 0])
basis(r⃗) = (normalize(r⃗), rotation*normalize(r⃗) )

include("atmosphere.jl")

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
	mean_anomaly,
	time_to_apoapsis
include("orbit.jl")


export Simulate!, SimulationParameters, standard_parameters
include("simulation.jl")


export plot_orbit, colTuple, animate_path
include("plotting.jl")

export Maneuver, runManeuver!, burn_time
include("maneuvers.jl")

export Model, runModel, reward
include("MachineLearning/modeling.jl")

#  going to have a naming conflict with `train_model` soon
export ∇reward, train_model! 
include("MachineLearning/gradient_descent.jl")

export train_model!, EpsilonExplorer
include("MachineLearning/reinforcement-learning.jl")

end # module AscentOptimization
