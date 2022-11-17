### A Pluto.jl notebook ###
# v0.19.14

using Markdown
using InteractiveUtils

# ╔═╡ 9a6fc386-5940-11ed-0eca-ef068c00402d
begin
	using Pkg;
	Pkg.activate("..");
	using StaticArrays
	using Plots
	using LinearAlgebra
	using Polynomials
	using Printf
	using GeometryBasics: Point2
	using AscentOptimization
	using LaTeXStrings
	using Base.Threads
	md"# Orbital Mechanics Sim"
end

# ╔═╡ 53d2281b-6283-4c7d-ab72-623af3044807
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

# ╔═╡ 6b93bee8-cd24-401c-aa59-98db5eeeb809
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

# ╔═╡ e7b64e5d-051e-4cb7-b076-acc4db6ea3ae
function Simulate!(ship::Ship, body::Planet; Δt)
	return if norm(ship.position) < body.radius
		#throw "you crashed"
		stick_to_ground!(ship, body, Δt)
	else
		Eulars!(ship, body; Δt)
	end
end

# ╔═╡ 6b8b106f-b704-4bf6-ba6b-9e53e1613426
function plot_orbit(ship, body; legend=false, aspect_ratio=1, kwargs...)
	
	ê = normalize(eccentricity_vector(ship,body))
	ϕ = acos(ê[1:2]⋅[1,0])
	θ = range(0,2π, length=100)
	R = body.radius
	#o = orbit(ship, body)
	apoapsis_point = Point2(ê[1:2] * apoapsis(ship,body))
	periapsis_point = Point2(-ê[1:2] * periapsis(ship,body))
	
	plt = plot(;aspect_ratio, legend, kwargs...)
	# plot!(plt, -o.(θ) .* cos.(θ  .- ϕ), -o.(θ) .* sin.(θ .- ϕ))
	plot!(plt, R .* cos.(θ), R .* sin.(θ))
	
	scatter!(plt, [apoapsis_point, periapsis_point])
	scatter!(plt, Point2(ship.position[1:2]))
	annotate!(plt, [
		(apoapsis_point...,"$(@sprintf("%.2f",apoapsis(ship,body)))"),
		(periapsis_point...,"$(@sprintf("%.2f",periapsis(ship,body)))")
	])
	return plt
end

# ╔═╡ 20707feb-a1b2-4050-83e3-b9ca5a1fd42f
let;
	struct EpsilonExplorer{A}
		ε::Float64  # not the specific energy
		learning_rate::Float64
		Q::Dict{A, Float64}
		N::Dict{A, Int}
	end
	function EpsilonExplorer(ε,α, action_space)
		return EpsilonExplorer(ε, α, Dict(action_space.=>0), Dict(action_space.=>0))
	end
end

# ╔═╡ b95ab490-48e3-457b-a6f2-fd689a49ec6f
function choice(e::EpsilonExplorer{A})::A where A
	return if rand() < e.ε
		rand([action for (action, value) in e.Q if isnan(value) || value==maximum(values(e.Q))])
	else
		rand(collect(keys(e.Q)))
	end
end	

# ╔═╡ d39d38e0-ceae-4bfe-80ed-5f3aaec3a7dd
function update!(learner::EpsilonExplorer{A}, action::A, reward) where A
	learner.Q[action] += learner.learning_rate * (reward - learner.Q[action])
	learner.N[action] += 1
end

# ╔═╡ 8b0ac269-2d0c-420a-9fe4-4d147f1c5885
struct SimulationState
	ship::Ship
	body::Planet
	results::Array{Ship,1}
	actions::Vector{Pair{Tuple{Float64, Float64},Float64}}
end

# ╔═╡ c4f8f244-5662-47e7-a152-855bd63e7dd6
isterminated(s::SimulationState) = s.ship.wet_mass <= 0

# ╔═╡ 9d5e9143-c0ae-4393-8137-39b0d48a45dc
function reward(state::SimulationState; target_altitude=12) 
	return reward(state.ship, state.body; target_altitude)
end	

# ╔═╡ 6f930d80-b47e-4cf1-94bf-d10e98afa951
function reward(ship::Ship,body::Planet; target_altitude=12)
	a = max(0, apoapsis(ship,body))
	p = max(0,periapsis(ship,body))
	ε = specific_energy(ship,body)
	return 4target_altitude^2 - 2(a - target_altitude)^2 - (p - target_altitude)^2
end

# ╔═╡ 9ea79efc-4e27-40ed-8d73-63c8fe399ec2
reward(a::Float64, p::Float64, t::Float64) = 4t^2 - 2(a - t)^2 - (p - t)^2

# ╔═╡ 79cd5514-89fc-4fea-b1d4-f194cff71a3c
function apply!(state::SimulationState, (pitch, throttle))
	(body, ship) = (state.body, state.ship)
	μ = body.gravitational_parameter

	state.ship.declination = (pitch + state.ship.declination ) % 2π
	state.ship.throttle = max(0, min(1, state.ship.throttle+throttle))
	Simulate!(state.ship, body, Δt=0.001)
	push!(state.results, deepcopy(state.ship))
	push!(state.actions, (pitch, throttle)=>reward(state))
	return state
end

# ╔═╡ eae30ca9-c3af-4704-9688-12307e242869
function runProblem(learner, state; steps=1000, target_altitude=12)
	for step in 1:steps
		Action = choice(learner)
		apply!(state, Action)
		Reward = reward(state;target_altitude)
		update!(learner, Action, Reward)
		isterminated(state) && (@info "you've been terminated"; break)
	end
end

# ╔═╡ 33634cbb-d6ab-4323-af2e-e9a04068ba38
action_space = Iterators.product(range(-π/6,π/6, length=100),range(0,1,length=10))

# ╔═╡ 43c6507e-8c39-4d2c-85f3-6059c6b1546c
const earth = Planet(1e7,1e3,0.1)

# ╔═╡ 99765703-7f0d-4901-9df4-3354f7153f29
colTuple(v::Vec) = ([v[1]], [v[2]], [v[3]])

# ╔═╡ 79a8278d-4211-4191-b9a0-8ea8ee4767ce
colTuple(v::Vector) = tuple(([u] for u in v)...)

# ╔═╡ 7f39761f-d438-461e-8c19-776b7cb4ec35
colTuple([100.0,100.0,100])

# ╔═╡ d1c4db1a-7564-4245-8658-4ede6b2e2364
weight(ship, body) = mass(ship) * surface_gravity(body)

# ╔═╡ 0b551d1f-7175-4042-b079-9d7047826964
twr(ship,body) = thrust(ship)/weight(ship,body)

# ╔═╡ a75031c3-5798-4b5e-9c8d-db44b6422ada
function heading(ship)
	(ρ̂, τ̂) = basis(ship.position)
	ϕ = ship.declination
	return ρ̂*cos(ϕ) + τ̂*sin(ϕ)
end

# ╔═╡ c60f30a4-91f7-4cce-b783-897d692f01e6
 body =  Planet(1e8,1e3,0.0)

# ╔═╡ 608bb9b6-0e15-4869-bdba-52abb51450d1
let g = surface_gravity(body)
global initial_ship =
	Ship(
		position=Vec(0, body.radius + 0.01,0),
		velocity=Vec(body.radius * body.angular_speed, 0.1,0),
		throttle=1,
		declination=0.0,
		dry_mass=20,
		max_thrust=2g*520,
		exhaust_velocity=350
	)
	initial_ship.wet_mass = initial_ship.max_thrust/2g - initial_ship.dry_mass
end

# ╔═╡ 759e3ed7-827e-4576-8a7a-3828a25b54c5
target_altitude = body.radius * 1.5

# ╔═╡ 29beddaf-acb1-4bf6-ba73-8cb1aaad3937
Base.@kwdef struct Maneuver
	stop_condition::Function
	throttle::Function
	declination::Function
end

# ╔═╡ 41a27fa4-ad7c-409d-80dc-cf431dea250c
function runManeuvers!(ship, maneuvers;  body, Δt=1e-6, step_limit=10^5)
	g = surface_gravity(body)	
	path = [deepcopy(ship)]
	Δt₀ = Δt 
	time_since_last_frame=0
	fps=20
	for (m,manuver) in enumerate(maneuvers)
		n=0
		while !manuver.stop_condition(ship)
			delta_v(ship) <=0 && throw("out of feul $m")
			n > step_limit && throw("iteration limit exceeded $m")
			specific_energy(ship, body) > 0 && throw("ship escaping $m")
			norm(ship.position) < body.radius && throw("ship crashed")
			ship.throttle = manuver.throttle(ship, Δt)
			ship.declination = manuver.declination(ship, Δt)
			Δt=min(Δt₀, 1/ norm(ship.velocity))
			Simulate!(ship,body; Δt)
			time_since_last_frame += Δt
			if time_since_last_frame > 1/fps
				push!(path, deepcopy(ship))
				time_since_last_frame = 0
			end
			n+=1
		end
	end
	@label done

	ship.throttle=0
	for t in 0:Δt:2
		Simulate!(ship,body;Δt)
		time_since_last_frame += Δt
		if time_since_last_frame > 1/fps
			push!(path, deepcopy(ship))
			time_since_last_frame = 0
		end
	end
	push!(path, deepcopy(ship))
	return path
	
end

# ╔═╡ 2c10ddd0-42bd-4938-b364-1e119f96c240
let N=10
	global Ω = [0] #range(0.445,0.4475,length=N)
	global deltaV=zeros(Float64,N) .- 10

	@threads for (i,ω) in collect(enumerate(Ω))
		maneuvers = [
		Maneuver( # initial burn
			stop_condition = (ship)->apoapsis(ship,body) >= target_altitude+0.1,
			throttle = (_,_) ->1,
			declination = (ship,Δt) ->  0 #min(π/2, ship.declination + ω*Δt)
		),
		Maneuver( # coast to apsis
			stop_condition = ship -> norm(ship.position) > target_altitude,
			throttle = (_,__) -> 0,
			declination = (_...) -> π/2
		),
		Maneuver( #rais apoapsis to appropriate height
			stop_condition = (ship)->apoapsis(ship,body) >= target_altitude,
			throttle = (_...) ->1,
			declination = (_...) -> π/2
		),
		Maneuver( #circularize
			stop_condition = (ship)->periapsis(ship,body) >= target_altitude,
			throttle = (_...) ->1,
			declination = (_...) -> π/2
		)
	]
		try
			(initial, final) = runManeuvers!(
				deepcopy(initial_ship),
				maneuvers;  
				body,
				Δt=1e-7,
				step_limit=1e8
			)[[1,end]]
			deltaV[i] = delta_v(final)
		catch  e
			@info "error" e ω
		end
	end
end

# ╔═╡ 4b1000cf-92dc-41ee-a755-acb01e043622
plot(Ω, max.(1, deltaV))

# ╔═╡ 6963d59b-a7f3-41a2-8bb2-51da22445a47
path = let	g = surface_gravity(body)
	ship=deepcopy(initial_ship)
	
	path = [deepcopy(ship)]
	
	ω = 0.44625
	Δt₀ = 1e-6
	Δt = Δt₀
	time_since_last_frame=0
	fps=10
maneuvers = [
		Maneuver( # initial burn
			stop_condition = (ship)->apoapsis(ship,body) >= target_altitude+0.1,
			throttle = (_,_) ->1,
			declination = (ship,Δt) ->  min(π/2, ship.declination + ω*Δt)
		),
		Maneuver( # coast to apsis
			stop_condition = ship -> norm(ship.position) > target_altitude,
			throttle = (_,__) -> 0,
			declination = (_...) -> π/2
		),
		Maneuver( #circularize
			stop_condition = (ship)->periapsis(ship,body) >= target_altitude,
			throttle = (_...) ->1,
			declination = (_...) -> π/2
		)
	]
	t=0
	for manuver in maneuvers
		while  !manuver.stop_condition(ship) && t < 10 
			delta_v(ship) <=0 && (@info "out of feul" ; @goto done)
			ship.throttle = manuver.throttle(ship,Δt)
			ship.declination = manuver.declination(ship,Δt)
			Δt=min(Δt₀, 1/ norm(ship.velocity))
			Simulate!(ship,body; Δt)
			time_since_last_frame += Δt
			t += Δt
			if time_since_last_frame > 1/fps
				push!(path, deepcopy(ship))
				time_since_last_frame = 0
			end
		end
	end
	@label done

	ship.throttle=0
	for t in 0:Δt:period(ship,body)
		Simulate!(ship,body;Δt)
		time_since_last_frame += Δt
		if time_since_last_frame > 1/fps
			push!(path, deepcopy(ship))
			time_since_last_frame = 0
		end
	end
	push!(path, deepcopy(ship))
	@info "" delta_v(ship)
	path
end

# ╔═╡ 9959a909-7fe1-48cf-b72c-98d38fa85da6
let ship = last(path)
	@info "" delta_v(ship) ship.velocity thrust_vector(ship)⋅ship.position apoapsis(ship, body) periapsis(ship,body) specific_energy(ship,body)

	plot(
		plot(reward.(path, Ref(body); target_altitude), label="reward"),
		plot(specific_energy.(path, Ref(body)), label="ε"),
		plot(delta_v.(path), label="Δv"),
		plot(twr.(path, Ref(body)), label=L"\frac{thrust}{weight}"),
		plot([ship.declination/π for ship in path], label="declination")
	)
end

# ╔═╡ 378e2614-3bad-4829-bf4a-55bc5dfe4161
function path_plot!(plt, path; kwargs...)
	points = Point2.(path[begin:end-1])
	δ = diff(path)
	U = [v[1] for v in δ]
	V = [v[2] for v in δ ]
	quiver!(plt, points; quiver=(U,V), kwargs...)
end

# ╔═╡ 4a7d25d3-26d3-4d0a-add0-773f86453068
function apsis_plot!(plt, path, a=0:1.2target_altitude, b = 0:1.2target_altitude)
	R = [reward(x,y, target_altitude) for x in a, y in b]
	
	contour!(plt, a,b,R,aspect_ratio=1)
	scatter!(plt, Point2(target_altitude,target_altitude),legend=false)
	subpath = path[1:100:end]
	path_plot!(
		plt, 
		collect.(AscentOptimization.apsies.(path; body))
	)
end

# ╔═╡ 3def2b5c-2838-4c8c-b237-2944361b3b7b
apsis_plot(args...) = apsis_plot!(plot(), args...)

# ╔═╡ 5ca79d0a-7c06-4f99-bc07-ff35f672a694
apsis_plot(path)

# ╔═╡ 48e3137b-180a-4524-becf-e80a3906aa53
let ship = deepcopy(initial_ship), g = surface_gravity(body)
	plt = plot()
	for r in range(0.01, body.radius*2, length=10)
		v_orbit = circular_orbit_speed(r, body)
		path = [Ship(
			position=Vec(0, r,0),
			velocity=v .* Vec(1, 0,0),
			throttle=1,
			declination=0.0,
			dry_mass=20,
			max_thrust=2g*520,
			exhaust_velocity=350
		) for v in range(0,v_orbit/1.2, length=10)]
		apsis_plot!(plt, path)
	end
	plt
end
# 43804 central station dr apartment 232 ashburn virginia 20147

# ╔═╡ 7de92c85-33ff-4369-aa83-26397c9e471f
let
anim = @animate for (n,ship) in enumerate(path)
		plt1 = plot_orbit(ship,body)
		plot!(Point2.(ship.position[1:2] for ship in path[1:n]))
		quiver!(Point2(ship.position[1:2]), quiver=colTuple(100*heading(ship)[1:2]))
		e =eccentricity_vector(ship,body)*100
		quiver!([0],[0], quiver=([e[1], [e[2]]]))
		
		plt2 = plot(Point2.(ship.velocity[1:2] for ship in path[1:n]))
		scatter!(Point2(last(path).velocity[1:2]))

		plt3 = plot(Point2.(thrust_vector(ship)[1:2] for ship in path[1:n]))
		scatter!(Point2(thrust_vector(last(path))[1:2]))
		
		plot(plt1,plt2, plt3)
end 
mov(anim)
end

# ╔═╡ Cell order:
# ╠═9a6fc386-5940-11ed-0eca-ef068c00402d
# ╠═53d2281b-6283-4c7d-ab72-623af3044807
# ╠═6b93bee8-cd24-401c-aa59-98db5eeeb809
# ╠═e7b64e5d-051e-4cb7-b076-acc4db6ea3ae
# ╠═6b8b106f-b704-4bf6-ba6b-9e53e1613426
# ╠═20707feb-a1b2-4050-83e3-b9ca5a1fd42f
# ╠═b95ab490-48e3-457b-a6f2-fd689a49ec6f
# ╠═d39d38e0-ceae-4bfe-80ed-5f3aaec3a7dd
# ╠═8b0ac269-2d0c-420a-9fe4-4d147f1c5885
# ╠═c4f8f244-5662-47e7-a152-855bd63e7dd6
# ╠═79cd5514-89fc-4fea-b1d4-f194cff71a3c
# ╠═9d5e9143-c0ae-4393-8137-39b0d48a45dc
# ╠═6f930d80-b47e-4cf1-94bf-d10e98afa951
# ╠═9ea79efc-4e27-40ed-8d73-63c8fe399ec2
# ╠═eae30ca9-c3af-4704-9688-12307e242869
# ╠═33634cbb-d6ab-4323-af2e-e9a04068ba38
# ╠═43c6507e-8c39-4d2c-85f3-6059c6b1546c
# ╠═e0672db1-5a3a-4336-b4bb-63030172974d
# ╠═99765703-7f0d-4901-9df4-3354f7153f29
# ╠═79a8278d-4211-4191-b9a0-8ea8ee4767ce
# ╠═7f39761f-d438-461e-8c19-776b7cb4ec35
# ╠═d1c4db1a-7564-4245-8658-4ede6b2e2364
# ╠═0b551d1f-7175-4042-b079-9d7047826964
# ╠═a75031c3-5798-4b5e-9c8d-db44b6422ada
# ╠═c60f30a4-91f7-4cce-b783-897d692f01e6
# ╠═608bb9b6-0e15-4869-bdba-52abb51450d1
# ╠═759e3ed7-827e-4576-8a7a-3828a25b54c5
# ╠═29beddaf-acb1-4bf6-ba73-8cb1aaad3937
# ╠═41a27fa4-ad7c-409d-80dc-cf431dea250c
# ╠═2c10ddd0-42bd-4938-b364-1e119f96c240
# ╠═4b1000cf-92dc-41ee-a755-acb01e043622
# ╠═6963d59b-a7f3-41a2-8bb2-51da22445a47
# ╟─9959a909-7fe1-48cf-b72c-98d38fa85da6
# ╟─378e2614-3bad-4829-bf4a-55bc5dfe4161
# ╟─4a7d25d3-26d3-4d0a-add0-773f86453068
# ╠═3def2b5c-2838-4c8c-b237-2944361b3b7b
# ╠═5ca79d0a-7c06-4f99-bc07-ff35f672a694
# ╟─48e3137b-180a-4524-becf-e80a3906aa53
# ╠═7de92c85-33ff-4369-aa83-26397c9e471f
