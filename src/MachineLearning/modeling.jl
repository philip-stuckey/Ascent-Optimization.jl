using Interpolations

"the model for the ascent uses a vector to generage a cubic spline of the angle of the ship from the verticle"
Base.@kwdef struct Model
	pitch_rate::Float64
	declination::Vector{Float64}
end

function pitch_times(model::Model)
	θ₀ = first(model.declination)
	θ₁ = last(model.declination)
	ω = abs(model.pitch_rate)
	t₀ = 0
	t₁ = (θ₁ ≈ θ₀) ? 1.0 : abs(θ₁ - θ₀) / ω
	
	return range(t₀, t₁, length=length(model.declination))
end	

"This function takes a mdoel and returns a new function that defines 
the declination (the angle of the ship from the vertical) over time.
If the time is outside the interval defined in the model, θₘₐₓ is returned"
declination(m::Model)  = extrapolate(
	scale(
		interpolate(
			m.declination, 
			BSpline(Cubic(Line(OnGrid())))
		),
		pitch_times(m)
	), last(m.declination)
)

"This function takes a model and returns a function that specfies the throttle of a
ship at a given time. right now it's clamped to full throttle"
throttle(::Model) = (_...)-> 1.0 #cubic_spline_interpolation(m.t, m.T)


# plugging the model into the maneuver system is relativly easy
function Maneuver(model::Model; done)
	tfunc = throttle(model)
	dfunc = declination(model)
	return Maneuver(
		done, 
		(_,t)->tfunc(t), 
		(_,t)->dfunc(t)
	)
end

"
	runModel(model::Model, ...)

take a model and some initial conditions, use the model to control the ascent of the
ship off the body, then coast to apopapsis and \"circularize\" (or at least get the periapsis over the `target_altitude`). 

The coast and circularize maneuvers requare their own meathods of `runManeuver` (see below)"
function runModel(model::Model, parameters::SimulationParameters; path=Ship[])
	
	ship = deepcopy(parameters.initial_ship)
	path === nothing || push!(path,deepcopy(ship))

	target_altitude = parameters.target_altitude
	margin = parameters.margin
	body = parameters.body

	ascent = Maneuver(
		model,
		done=(_...)->apoapsis(ship, body) >= target_altitude + margin,
	)

	coast = Maneuver(
		done=(_...)->norm(ship.position) >= target_altitude,
		throttle=0.0,
		declination=π/2
	)

	circularize = Maneuver(
		done=(_...)->periapsis(ship, body) >= target_altitude,
		throttle=1.0,
		declination= π/2
	)

	for maneuver in (ascent, coast, circularize)
		runManeuver!(ship, maneuver, parameters; path)
	end

	return (ship, path)
end

function reward(θ, parameters=standard_parameters)
	model =Model(θ[1], θ[2:end])
	return reward(model, parameters)
end

function reward(model::Model, parameters=standard_parameters)
	body = parameters.body
	target_altitude = parameters.target_altitude
	
	(ship, _) = runModel(model, parameters, path=nothing)
	return delta_v(ship)*(target_altitude < periapsis(ship,body))
end
