function Eulars!(ship::Ship, body::Planet; Δt=0.001)
	μ = body.gravitational_parameter
	(ρ̂, _) = basis(ship.position)
	r⃗ = ship.position
	m = mass(ship)
	v = ship.velocity
	ϕ = ship.declination
	T⃗ = thrust_vector(ship)
	ṁ = mass_flow_rate(ship)
	
	ship.wet_mass -= ṁ*Δt
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
	T⃗ -= min(0,T⃗⋅ρ̂)*ρ̂  # remove the downward component ot T⃗
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