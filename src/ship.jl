using Unitful: kg, m, s, kN

Base.@kwdef mutable struct Ship
    const dry_mass::Mass
	const exhaust_velocity::Velocity = 1000m/s
	const max_thrust::Force = 100kN
		
	throttle::Float64=1
	declination::Float64 = 0
	
	fuel_mass::Mass = 0kg  # we're out of cream, would you like "no milk" instead?
    position::Vec{L}  where L <:Length
    velocity::Vec{T} where T <:Velocity
end

function Ship(
	body::Planet;
	Δv::Velocity, 
	dry_mass::Mass, 
	TWR=2, 
	vₑₓ=1000m/s,
	position=Vec(0m,body.radius,0m),
	velocity=Vec(surface_speed(body),0m/s,0m/s)
)
	g = surface_gravity(body)
	total_mass = dry_mass*exp(Δv/vₑₓ)
	return Ship(;
		dry_mass,
		exhaust_velocity=vₑₓ,
		max_thrust = TWR*g*total_mass,
		fuel_mass = total_mass - dry_mass,
		position,
		velocity
	)
end

function Base.setproperty!(ship::Ship, name::Symbol, x) 
	v = if name == :throttle
		clamp(x, zero(x), oneunit(x))
	elseif name == :feul_mass
		max(zero(x), x)
	else
		x
	end
	setfield!(ship, name, v)
end

mass(ship::Ship) = ship.dry_mass + ship.fuel_mass
thrust(ship::Ship) = ship.max_thrust * ship.throttle * (ship.fuel_mass > 0kg)

function heading(ship)
	(ρ̂, τ̂) = basis(ship.position)
	ϕ = ship.declination
	return ρ̂*cos(ϕ) + τ̂*sin(ϕ)
end

function thrust_vector(ship)::Vec{Force} 
	return thrust(ship) * heading(ship)
end

mass_flow_rate(ship) = thrust(ship) / ship.exhaust_velocity

function delta_v(ship)
	vₑₓ = ship.exhaust_velocity
	m₀ = ship.fuel_mass
	m₁ = ship.dry_mass
	return vₑₓ * log(1 + m₀/m₁)
end

specific_angular_momentum(ship) = (ship.position × ship.velocity)