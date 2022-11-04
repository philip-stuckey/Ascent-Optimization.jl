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

function period(ship, body)
	a = semi_major_axis(ship,body)
	μ = body.gravitational_parameter
	return 2π*√(a^3 / μ)
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
    return  if ε < 0
		(-μ ± √(μ^2 + 2ε*(l⋅l))) ./ 2ε
	else
		((-μ + √(μ^2 + 2ε*(l⋅l))) ./ 2ε, Inf)
	end
end

apoapsis(ship,body) = apsies(ship;body)[2]
periapsis(ship,body) = apsies(ship;body)[1]

function circular_orbit_speed(ship::Ship, body::Planet) 
	r = norm(ship.position)
	μ = body.gravitational_parameter
	return √(μ/r)
end


function orbit(ship, body)
	l⃗ = specific_angular_momentum(ship)
	e = eccentricity(ship, body)
	μ = body.gravitational_parameter
	return θ -> (l⃗⋅l⃗) / (μ * (1 + e * cos(θ)))
end

function orbital_angular_speed(ship,body)
	T = period(ship,body)
	return 2π/T
end

function time_to_apsis(ship,body)
	T= period(ship,body)
	θ = ship.position ⋅ eccentricity_vector(ship,body)
	return θ*T/(2π)
end

weight(ship, body) = mass(ship) * surface_gravity(body)
twr(ship,body) = thrust(ship)/weight(ship,body)