### A Pluto.jl notebook ###
# v0.19.14

using Markdown
using InteractiveUtils

# ╔═╡ e3498a58-5f72-11ed-29ea-6fea4710b976
begin
	using Pkg
	Pkg.activate(".")
	#include("src/AscentOptimization.jl")
end

# ╔═╡ 742b9989-7121-4e99-b518-3cdad6d3349d
begin
	using AscentOptimization
	using LinearAlgebra
end

# ╔═╡ 68ace98c-49d5-43e2-8d5b-231232a4c111
using Plots, Printf

# ╔═╡ 6a051cea-e9c6-43e7-a9dc-6b84dc94cf02


# ╔═╡ d64ce408-fdc0-44f6-8095-a924029a7b22
runModel(Model(0.4, [0.0, π/2]), standard_parameters)

# ╔═╡ 26cd5fee-fc67-47d3-a73f-d4ee9c88e882
let target_altitude=15000
	params = SimulationParameters(
		initial_ship = standard_parameters.initial_ship,
		body = standard_parameters.body,
		margin = standard_parameters.margin,
		target_altitude = standard_parameters.target_altitude,
		time_step = standard_parameters.time_step,
		time_limit=13.0625
	)
	body = params.body
	target_altitude = params.target_altitude
	ship = deepcopy(params.initial_ship)
	
	ascent = Maneuver(
		Model(0.2, [0.0, π/2]),
		done=(_...)->apoapsis(ship, body) >= target_altitude,
	)
	global path=Ship[]
	runManeuver!(ship, ascent, standard_parameters; path)
	circularize =  Maneuver(
		(ship1)->periapsis(ship1, body) >= target_altitude,
		(ship,_...)->float(time_to_apoapsis(ship,body) <= 0.5 ),
		(ship,_...)->declination=π/2 - 0.05* (norm(ship.position) < target_altitude)
	)
	runManeuver!(ship, circularize, params; path)
end

# ╔═╡ ff1d8c9e-c135-477e-84d8-142f23cbe35c
animate_path(path, standard_parameters.body)

# ╔═╡ 0c431645-c1e7-44ee-95a6-8baa9c809198
plot(eccentricity.(path, Ref(standard_parameters.body)))

# ╔═╡ bcd4f0fb-4af6-4d1c-a0cd-80f4932e95d2
let
	plt1 = plot(
		time_to_apoapsis.(path, Ref(standard_parameters.body)), 
		label="time to apoapsis"
	)
	plt2 = plot(norm.(ship.position for ship in path))
	plot(plt1, plt2, layout=(2,1))
end

# ╔═╡ 86075f02-9eed-4e36-a73c-99156ad5419f
let body = standard_parameters.body
	ship = Ship(\Delta)

# ╔═╡ fe094aee-53f1-40e3-aac4-7779415d6a99
apoapsis(ship, standard_parameters.body)

# ╔═╡ 6bdc6bc9-abf1-4b88-9dba-7942735fac34
periapsis(ship, standard_parameters.body)

# ╔═╡ d79f2afd-341d-49e5-a9b0-eebce055ed58
delta_v(ship)

# ╔═╡ fd65d988-22c5-4953-80a1-359a97b4712a
gr()

# ╔═╡ 384ee66f-2c34-48dc-90b1-1954286e1b46
let body = standard_parameters.body
	#path=path[1:200]
	plot(
		#yformatter=x->@sprintf("%0.2fπ", x/π),
		legend=:bottomright
	)
	#plot!(true_anomaly.(path, Ref(body)), label="true anomaly")
	#plot!(eccentric_anomaly.(path, Ref(body)),label="eccentric anomaly")	
	#plot!(mean_anomaly.(path, Ref(body)), label="mean anomaly")
	
	plot!(time_since_periapsis.(path, Ref(body)) , label="time since periapsis")
	plot!(time_to_apoapsis.(path, Ref(body)), label="time to apoapsis")
	#plot!(period.(path, Ref(body)), label="period")
	
	#plot!(eccentricity.(path, Ref(body)), label="eccentricity")
	#plot!(cos.(true_anomaly.(path, Ref(body))), label="cos(true anomaly)")

	(_, n) = findmax(norm.(ship.position for ship in path))
	vline!([n], label="apoapsis reached")
	hline!([1.7011644608776602], label=nothing)
	#plot!([ship.throttle for ship in path], label="throttle")
end

# ╔═╡ c69ba80c-31d8-4258-81d3-f09e3402a5e4
# ╠═╡ disabled = true
#=╠═╡
animate_path(filter!(s->norm(s.position)>1000,path), standard_parameters.body)
  ╠═╡ =#

# ╔═╡ Cell order:
# ╠═e3498a58-5f72-11ed-29ea-6fea4710b976
# ╠═742b9989-7121-4e99-b518-3cdad6d3349d
# ╠═6a051cea-e9c6-43e7-a9dc-6b84dc94cf02
# ╠═d64ce408-fdc0-44f6-8095-a924029a7b22
# ╠═26cd5fee-fc67-47d3-a73f-d4ee9c88e882
# ╠═ff1d8c9e-c135-477e-84d8-142f23cbe35c
# ╠═0c431645-c1e7-44ee-95a6-8baa9c809198
# ╠═bcd4f0fb-4af6-4d1c-a0cd-80f4932e95d2
# ╠═86075f02-9eed-4e36-a73c-99156ad5419f
# ╠═fe094aee-53f1-40e3-aac4-7779415d6a99
# ╠═6bdc6bc9-abf1-4b88-9dba-7942735fac34
# ╠═d79f2afd-341d-49e5-a9b0-eebce055ed58
# ╠═68ace98c-49d5-43e2-8d5b-231232a4c111
# ╠═fd65d988-22c5-4953-80a1-359a97b4712a
# ╠═384ee66f-2c34-48dc-90b1-1954286e1b46
# ╠═c69ba80c-31d8-4258-81d3-f09e3402a5e4
