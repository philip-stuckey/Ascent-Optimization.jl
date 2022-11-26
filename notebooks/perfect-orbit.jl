### A Pluto.jl notebook ###
# v0.19.15

using Markdown
using InteractiveUtils

# ╔═╡ 313aad86-6d1b-11ed-072e-430b9066a8b0
begin
	import Pkg
	Pkg.activate("..")
end

# ╔═╡ b63b1c33-0358-461c-97c3-fe0c6601d182
using AscentOptimization

# ╔═╡ c1f0a09f-369d-4c23-b2ab-36229a30de42
using Unitful: kg, m, s

# ╔═╡ 78a9d6b8-c394-43ab-a4d7-1379729c165d
using LinearAlgebra

# ╔═╡ 5ddc823e-474e-489e-8cd0-1e9c6021c78d
using Plots

# ╔═╡ 024ad44a-32fd-418f-8d4d-c05352b903ee
begin
	import Unitful
	using Unitful: Quantity
	const 𝕃 = Unitful.𝐋
	const 𝕋 = Unitful.𝐓
end

# ╔═╡ f0326ea8-18d0-4803-be71-69bbba5761c1
include("../src/planets.jl")

# ╔═╡ 1cd88cab-5c11-4b9b-bf30-d5559fa6b71d
params = SimulationParameters(
	body=Gilly,
	initial_ship = Ship(
		Gilly, 
		position=Vec(0m, Gilly.radius, 0m),
		Δv=2000.0m/s,
		dry_mass=10.0kg, 
		TWR=3
	),
	margin = 1m,
	target_altitude = Gilly.radius*1.2,
	snapshot_rate=1/(1s),
	time_limit=3600s,
	time_step=1e-3s
)

# ╔═╡ ac816a6c-0d3c-4129-9a39-35f803f3383a
weight(ship,body)=mass(ship)*body.gravitational_parameter/norm(ship.position)^2

# ╔═╡ fc5d8ead-6a3f-4237-896b-6984223a7c41
function basic_ascent!(
	params=params,
	ascent_maneuver= AscentOptimization.Maneuver(
		(ship, _...)->AscentOptimization.apoapsis(ship, params.body) >= params.target_altitude+params.margin,
		(_...)-> 1.0,
		(ship, _...) ->atan(thrust(ship), weight(ship, params.body)+0.25kg*m/s^2)
	), 
	coast_maneuver= AscentOptimization.Maneuver(
		done=(ship, _...)->norm(ship.position) >= params.target_altitude,
		throttle=0.0,
		declination=π/2
	), 
	circularize_maneuver = AscentOptimization.Maneuver(
		done=(ship, _...)->AscentOptimization.periapsis(ship, params.body) >= params.target_altitude,
		throttle=1.0,
		declination= π/2
	);
	path=nothing
)

	ship = deepcopy(params.initial_ship)
	for manuever in (ascent_maneuver, coast_maneuver, circularize_maneuver)
		runManeuver!(ship, manuever, params; path)
	end
	return (ship,path)
end

# ╔═╡ 93b38985-7287-4915-b307-85eef8947b22
let ship = params.initial_ship, body=params.body
	atan(thrust(ship), weight(ship, body))
end

# ╔═╡ f77ac3d7-d98a-4b35-8dca-d3992cc5154d
let (ship, path) = basic_ascent!(params, path=Ship[])
	plot(diff(norm.(ship.position for ship in path)) .> zero(AscentOptimization.Length))
	#xlims!((0,100))
end

# ╔═╡ 3010c1e4-4fc8-42e0-90eb-c63d4ddd2851
(ship, path) = basic_ascent!(params, path=Ship[])

# ╔═╡ a388e75a-9401-407d-9330-30ca9ebe95e5
animate_path(
		path[begin:5:end], 
		params.body, 
		#xlims=(-100m, 100m),
		#ylims=(13000m, 13100m)
	)

# ╔═╡ ef5bbbe0-6765-486b-bbeb-a6b651b747f6
abstract type AbstractScalarOrbit end

# ╔═╡ f48ec539-a224-4c80-9b26-7838b7d08507
struct ScalarOrbit <: AbstractScalarOrbit
	gravitational_parameter::Quantity{Float64, 𝕃^3*𝕋^-2, typeof(m^3/s^2)}
	apoapsis::AscentOptimization.Length
	periapsis::AscentOptimization.Length
	ScalarOrbit(μ, apsis1, apsis2)=new(μ, max(apsis1,apsis2), min(apsis1,apsis2))
end

# ╔═╡ a98db8de-9d23-4812-aa0a-8b0f8fcf3c4c
gravitational_parameter(o::ScalarOrbit)=o.gravitational_parameter

# ╔═╡ 19076259-d82c-43d0-a36b-613891c724fe
periapsis(o::ScalarOrbit)=o.periapsis

# ╔═╡ 14b2f33d-be88-4ac7-99e8-77a4365a7fa7
apoapsis(o::ScalarOrbit)=o.apoapsis

# ╔═╡ ea35b9dc-0ff8-47c7-a0f8-af1ac9f02951
struct CircularOrbit <: AbstractScalarOrbit
	gravitational_parameter::Quantity{Float64, 𝕃^3*𝕋^-2, typeof(m^3/s^2)}
	radius::AscentOptimization.Length
end

# ╔═╡ 305c5046-b093-4067-a1bb-0c128021d2bf
apoapsis(o::CircularOrbit) = o.radius

# ╔═╡ bf96bcbe-57e2-46f3-bb80-def96993c47b
periapsis(o::CircularOrbit) = o.radius

# ╔═╡ 278c2d1f-e9ba-4e0d-9b80-3ac29eedde95
function periapsis_velocity(o::AbstractScalarOrbit)
	μ = gravitational_parameter(o)
	ra = apoapsis(o)
	rp = periapsis(o)
	return sqrt(2μ*ra*(ra - rp)/(rp * (ra^2 - rp^2)))
end

# ╔═╡ 3abab020-8565-4eb6-b732-db944c6da740
function apoapsis_velocity(o::AbstractScalarOrbit)
	μ = gravitational_parameter(o)
	ra = apoapsis(o)
	rp = periapsis(o)
	return sqrt(2μ*rp*(ra - rp)/(ra * (ra^2 - rp^2)))
end

# ╔═╡ e20396d1-8406-46c6-8f42-ac471860a1d6
function apoapsis_velocity(o::CircularOrbit) 
	r = o.radius
	μ = o.gravitational_parameter
	return √(μ/r)
end

# ╔═╡ 67d72245-d67f-4115-8695-5dba3bc8bc14
periapsis_velocity(o::CircularOrbit) = apoapsis_velocity(o)

# ╔═╡ 74cc9ff9-76ce-4eb6-9f24-c0b03c318185
ideal_Δv = let body=params.body, r=params.target_altitude+params.margin
	μ = body.gravitational_parameter
	initial_orbit = ScalarOrbit(μ, 0m, body.radius)
	final_orbit = CircularOrbit(μ, r)
	transfer_orbit = ScalarOrbit(μ, body.radius, r)
	ascent_Δv = periapsis_velocity(transfer_orbit) - apoapsis_velocity(initial_orbit)
	circularize_Δv = periapsis_velocity(final_orbit)-apoapsis_velocity(transfer_orbit)
	delta_v(params.initial_ship) - ascent_Δv - circularize_Δv
end

# ╔═╡ 259f101c-72e6-495e-b2b7-42fa446211a3
md"""
ideal_Δv - delta_v(ship) = $(ideal_Δv - delta_v(ship))

**close enough**
"""

# ╔═╡ Cell order:
# ╠═313aad86-6d1b-11ed-072e-430b9066a8b0
# ╠═b63b1c33-0358-461c-97c3-fe0c6601d182
# ╠═f0326ea8-18d0-4803-be71-69bbba5761c1
# ╠═c1f0a09f-369d-4c23-b2ab-36229a30de42
# ╠═1cd88cab-5c11-4b9b-bf30-d5559fa6b71d
# ╠═78a9d6b8-c394-43ab-a4d7-1379729c165d
# ╠═ac816a6c-0d3c-4129-9a39-35f803f3383a
# ╠═fc5d8ead-6a3f-4237-896b-6984223a7c41
# ╠═93b38985-7287-4915-b307-85eef8947b22
# ╠═5ddc823e-474e-489e-8cd0-1e9c6021c78d
# ╠═f77ac3d7-d98a-4b35-8dca-d3992cc5154d
# ╠═3010c1e4-4fc8-42e0-90eb-c63d4ddd2851
# ╠═a388e75a-9401-407d-9330-30ca9ebe95e5
# ╠═74cc9ff9-76ce-4eb6-9f24-c0b03c318185
# ╠═259f101c-72e6-495e-b2b7-42fa446211a3
# ╠═ef5bbbe0-6765-486b-bbeb-a6b651b747f6
# ╠═024ad44a-32fd-418f-8d4d-c05352b903ee
# ╠═f48ec539-a224-4c80-9b26-7838b7d08507
# ╠═a98db8de-9d23-4812-aa0a-8b0f8fcf3c4c
# ╠═19076259-d82c-43d0-a36b-613891c724fe
# ╠═14b2f33d-be88-4ac7-99e8-77a4365a7fa7
# ╠═278c2d1f-e9ba-4e0d-9b80-3ac29eedde95
# ╠═3abab020-8565-4eb6-b732-db944c6da740
# ╠═ea35b9dc-0ff8-47c7-a0f8-af1ac9f02951
# ╠═305c5046-b093-4067-a1bb-0c128021d2bf
# ╠═bf96bcbe-57e2-46f3-bb80-def96993c47b
# ╠═e20396d1-8406-46c6-8f42-ac471860a1d6
# ╠═67d72245-d67f-4115-8695-5dba3bc8bc14
