
Base.@kwdef mutable struct Ship
    const dry_mass::Float64
	const exhaust_velocity::Float64 =1000
	const max_thrust::Float64 = 100
		
	throttle::Float64=1
	declination::Float64 = 0
	
	wet_mass::Float64 = 0
    position::Vec
    velocity::Vec
end

function Ship(
	body::Planet;
	Δv, 
	dry_mass, 
	TWR=2, 
	vₑₓ=1000,
	position=Vec(0,body.radius,0),
	velocity=Vec(surface_speed(body),0,0)
)
	g = surface_gravity(body)
	total_mass = dry_mass*exp(Δv/vₑₓ)
	return Ship(;
		dry_mass,
		exhaust_velocity=vₑₓ,
		max_thrust = TWR*g*total_mass,
		wet_mass = total_mass - dry_mass,
		position,
		velocity
	)
end

function Base.setproperty!(ship::Ship, name::Symbol, x) 
	v = if name == :throttle
		min(oneunit(x), max(zero(x), x))
	elseif name == :wet_mass
		max(zero(x), x)
	else
		x
	end
	setfield!(ship, name, v)
end

mass(ship) = ship.dry_mass + ship.wet_mass
thrust(ship) = ship.max_thrust * ship.throttle * (ship.wet_mass > 0)

function heading(ship)
	(ρ̂, τ̂) = basis(ship.position)
	ϕ = ship.declination
	return ρ̂*cos(ϕ) + τ̂*sin(ϕ)
end

function thrust_vector(ship) 
	return thrust(ship) * heading(ship)
end

mass_flow_rate(ship) = thrust(ship) / ship.exhaust_velocity

function delta_v(ship)
	vₑₓ = ship.exhaust_velocity
	m₀ = mass(ship)
	m₁ = ship.dry_mass
	return vₑₓ * log(m₀/m₁)
end

specific_angular_momentum(ship) = (ship.position × ship.velocity)