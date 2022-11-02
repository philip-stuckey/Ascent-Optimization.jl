module OrbitalMechanics

using StaticArrays

const Vec = SVector{3, Float64}
Base.@kwdef mutable struct Ship
    const dry_mass::Float64
	const exhaust_velocity::Float64 =1000
	const max_thrust::Int = 100
		
	throttle::Float64=0
	wet_mass::Float64 = 0
	declination::Float64 = 0
    position::Vec
    velocity::Vec
end

function Base.setproperty!(ship::Ship, name::Symbol, x) 
	if name == :throttle
		setfield!(ship,:throttle, min(1.0, max(0.0, x)))
	else
		setfield!(ship,name, x)
	end
end

mass(ship) = ship.dry_mass + ship.wet_mass
thrust(ship) = ship.max_thrust * ship.throttle * (ship.wet_mass > 0)



end # module OrbitalMechanics
