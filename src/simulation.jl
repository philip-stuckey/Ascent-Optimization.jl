using  Unitful: m,kg,s
Base.@kwdef struct SimulationParameters
    body::Planet
    initial_ship::Ship
    time_step::Time
	margin::Length
	target_altitude::Length
	snapshot_rate::Int=60
	time_limit::Time=20s
end


const default_body = Planet(1e8m^3/s^2,1e3m,0.0/s)

const standard_parameters = SimulationParameters(
	body=default_body,
	initial_ship=Ship(default_body,Δv=2000.0m/s, TWR=2, dry_mass=20.0kg),
	time_step=(10^-5)s,
	margin=10.0m,
	target_altitude=default_body.radius*1.5
)

function Eulars!(ship::Ship, body::Planet; Δt=0.001s)
	μ = body.gravitational_parameter
	(ρ̂, _) = basis(ship.position)
	r⃗ = ship.position
	m = mass(ship)
	v = ship.velocity
	T⃗ = thrust_vector(ship)
	ṁ = mass_flow_rate(ship)
	
	ship.fuel_mass -= ṁ*Δt
	a⃗ = T⃗/m - ρ̂ * μ/(r⃗⋅r⃗) - v*ṁ/m
	ship.velocity += a⃗*Δt
	ship.position += ship.velocity * Δt
	
	return ship
end


function stick_to_ground!(ship::Ship, body::Planet, Δt)
	ω = body.angular_speed
	R = body.radius
	
	(ρ̂, τ̂) = basis(ship.position)
	T⃗ = thrust_vector(ship)
	m = mass(ship)
	T⃗ -= min(zero(T⃗⋅ρ̂),T⃗⋅ρ̂)*ρ̂  # remove the downward component ot T⃗
	ṁ = mass_flow_rate(ship)
	
	ship.fuel_mass -= ṁ*Δt
	a⃗ = T⃗/m
	
	ship.velocity = ω*R * τ̂ + a⃗*Δt
	ship.position = normalize(ship.position)*R
	ship.position += ship.velocity*Δt
	return ship
end

function Simulate!(ship::Ship, body::Planet; Δt)
	return if norm(ship.position) < body.radius
		#throw "you crashed"
		stick_to_ground!(ship, body, Δt)
	else
		Eulars!(ship, body; Δt)
	end
end
