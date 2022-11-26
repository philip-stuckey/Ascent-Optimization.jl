### A Pluto.jl notebook ###
# v0.19.15

using Markdown
using InteractiveUtils

# ╔═╡ 6ddda634-6b63-11ed-0766-a57e14d8b741
begin
	using Pkg
	Pkg.activate("..")
end

# ╔═╡ b051ae45-acc7-4de8-8ad6-abcff4fe6e17
using AscentOptimization

# ╔═╡ 102329cb-8c56-43f3-a1a4-7fab9984576a
using Unitful: kg, m, s

# ╔═╡ ee3eb72f-741e-4e84-9809-346f30e9f84e
using ThreadsX

# ╔═╡ c88cbe4e-a4b6-4174-8d55-69efc4e8beaa
using LinearAlgebra

# ╔═╡ d11ab9f0-e8f9-46ee-85e4-596ff624cc65
using Plots

# ╔═╡ 84dc7ae7-0f7c-45d4-9541-813f93348797
using FiniteDiff: finite_difference_gradient

# ╔═╡ 537e3669-2861-4ac2-925c-a84a6f2d473a
using Statistics

# ╔═╡ e2527ad1-1742-478e-bc52-7d9c20e535af
include("../src/planets.jl")

# ╔═╡ 6c330af0-b2ea-4d8c-888f-b76356b38119
params = SimulationParameters(
	body=Gilly,
	initial_ship = Ship(
		Gilly, 
		position=Vec(0m, Gilly.radius+10m, 0m),
		Δv=2000.0m/s,
		dry_mass=10.0kg, 
		TWR=3
	),
	margin = 100m,
	target_altitude = Gilly.radius*1.2,
	snapshot_rate=1/(1s),
	time_limit=3600s,
	time_step=1e-3s
)

# ╔═╡ 1f0c5d25-cab5-4765-b0d5-473dc601778f
surface_speed(Gilly)

# ╔═╡ c2290d22-e4eb-42c6-b01b-c13287caaab0
let (ship,path) = runModel(Model(0.1/s, [1.39626,1.22173]), params)
	@info "" delta_v(ship)
	animate_path(path[begin:10:end], params.body)
end

# ╔═╡ 08226305-1421-4a82-bb80-5df629f59957
N=10

# ╔═╡ 700255bf-88f8-4b72-9168-7bcfebf11184
Ω = range(0.1, 3, length=N)

# ╔═╡ 615bab31-67a4-4f56-923f-4011bafcba18
Θ = range(0, π/2, length=N)

# ╔═╡ 8ba4a032-a56e-4d5c-8adf-fac41bf31d2c
rewards = ThreadsX.map(x-> reward(Vec(x), params), Iterators.product(Ω, Θ, Θ));

# ╔═╡ 30b8ca8e-d468-4afd-b0b7-f58b94b37246
∇rewards = ThreadsX.map(x-> norm(∇reward(collect(x), params)), Iterators.product(Ω, Θ, Θ));

# ╔═╡ 04e5b976-58d5-4cfc-b95c-6bc52bbb1cfc
finite_difference_gradient(
	x->reward(x, params), 
	[0.5, 0, 0.166π], 
	Base.Val{:central}, 
	AscentOptimization.RewardType; 
	relstep=0.01
)

# ╔═╡ 89f65794-4041-40a1-a860-dedbef7fdf23
let M = 8
	heatmap(Ω, Θ, rewards[:,:,M], title=string(Θ[M]/π)*"π")
end

# ╔═╡ 1854c288-9503-4d99-82f6-b316db3070ef
findmax(rewards)

# ╔═╡ 128a3180-5eee-4ffc-a53b-83833bb71ccc
Ω[1], Θ[9], Θ[8]

# ╔═╡ be0050c5-21eb-41ce-a2ea-18aaaf1d64a6
let M = 8
	heatmap(Θ, Θ,∇rewards[M,:,:], title=string(Ω[M]/π)*"π")
end

# ╔═╡ c8be7b9c-4000-4f9f-a654-eace228b77cd
findmax(∇rewards)

# ╔═╡ 58e0eb05-fe1b-4e59-bda6-b6f9bc95967d
(Ω[3], Θ[3], Θ[9])

# ╔═╡ 3f5a66b7-6b6d-4875-83be-8493e1d3724b
Θ[6]/π

# ╔═╡ 03af5d6b-826a-45d9-8396-50a84e6bc07f
reinfocement_model = let
	global reinforcement_rewards = AscentOptimization.RewardType[]
	train_model!(
		Model(0.3/s, [0, 2π]), 
		params, 
		EpsilonExplorer(1.0, 0.1); 
		rewards=reinforcement_rewards,
		steps=200
	)
end

# ╔═╡ 856364cf-4804-4c5a-b2d3-b0c72697d90d
αs = range(0,1, length=10)

# ╔═╡ d4f7bf69-ddb0-45d2-a895-6fb0794c5198
εs = range(0,1, length=10)

# ╔═╡ bb53229e-2a85-4762-a1fb-980897d105d7
reinforcement_models = ThreadsX.map(
	function (x)
		local rewards=AscentOptimization.RewardType[]
		local explorer=EpsilonExplorer(x[1],x[2])
		local model = train_model!(
			Model(0.3/s, [0, 2π]), 
			params, 
			explorer; 
			rewards=rewards,
			steps=350
		)
		return model, rewards, explorer
	end
	,
	Iterators.product(αs, εs)
)

# ╔═╡ 71876131-4153-4fb9-97be-145922cadd4f
findmax(x->last(x[2]),reinforcement_models)

# ╔═╡ de4659fc-2de4-44b6-bad9-4df78938f8d7
delta_v(params.initial_ship) - circular_orbit_speed(params.target_altitude, params.body)

# ╔═╡ 0245f903-9feb-4b5f-abf3-3f34bb815dd4
function transfer_velocity(r₀, r; body = params.body)
	μ = body.gravitational_parameter	
	return √(2μ*r*(r-r₀)/ (r₀*(r^2 - r₀^2)))
end

# ╔═╡ 01905a71-be04-4214-a226-4d4a4686a1d2
function transfer_deltaV(r₀, v₀, r; body::Planet = params.body)
	μ = body.gravitational_parameter
	vₜ = transfer_velocity(r₀, r; body)
	burn1_Δv = vₜ - v₀
	
	local apoapsis_velocity = √(vₜ^2 + 2μ*(1/r - 1/r₀))
	burn2_Δv = circular_orbit_speed(r, body) - apoapsis_velocity
	return burn2_Δv + burn1_Δv
end

# ╔═╡ d5628d11-8301-4844-833d-3d77e4ca69dd
let r₀ = 1000m, body=params.body, μ = body.gravitational_parameter, v₀=√(μ/r₀)
	rs = range(0m, 10000m, step=1m)
	escape_velocity=√(2μ/r₀)
	plot(rs, [transfer_deltaV(r₀,v₀, r; body) for r in rs])
	plot!(rs, [√(2μ/r) for r in rs])
	hline!([√(μ/r₀), √(2μ/r₀)])
	vline!([r₀])
	ylims!((-escape_velocity,escape_velocity))
end

# ╔═╡ 0330467c-59a2-4941-b881-42fb0893e85b
let ship=params.initial_ship, r₀ = norm(ship.position), v₀ = norm(ship.velocity)
	delta_v(ship) - transfer_deltaV(r₀, v₀, params.target_altitude; body=Gilly)
end

# ╔═╡ 73ee0b41-8368-4f1f-9806-ef97780cbe7b
begin
	plt =plot(
		#legend=:bottomright,
		legend=false
		#xlims=(250,300),
		#ylims=(1960,1970)
	)

	ship=params.initial_ship
	r₀ = norm(ship.position)
	v₀ = norm(ship.velocity)
	r = params.target_altitude
	body = params.body
	
	optimal_deltaV = delta_v(ship) - transfer_deltaV(r₀, v₀, r; body)
	for (model, rewards, explorer) in reinforcement_models
		if last(rewards) > 1966m/s
			plot!(
				plt,
				rewards ./ optimal_deltaV,
				#label=string(model)
				label=string(explorer)
			)
		end
	end
	vline!(plt, [375])
	plt
end

# ╔═╡ 2f621adb-7207-4831-8123-c081d66588a9
let ship=params.initial_ship, body=params.body
	Δv = delta_v(ship)
	twr = AscentOptimization.twr(ship, body)
	Δv * √(twr^2-1)/twr
end

# ╔═╡ 3947e65c-ebd2-49ed-b842-4a2343457620
last(reinforcement_models[1][2]) > (1965m/s)

# ╔═╡ 5dd533f3-b6c4-4591-ac3d-ec3960c6ece7
histogram([last(M[2]) for M in reinforcement_models[:]])

# ╔═╡ efc466a6-5caf-4fdd-836f-c93c168b88ae
let models = reinforcement_models
	heatmap(
		αs,
		εs,
		[maximum(M[2] .* (1-1e-5).^(1:length(M[2]))) for M in permutedims(models)]
	)
end

# ╔═╡ c87cbaf2-4061-4991-9263-737509f2a945
let plt = plot(legend=false, xlims=(0,400))
	for model in reinforcement_models
		plot!(plt, [M for M in model[2]] .* (1-1e-6) .^(1:length(model[2])))
	end	
	plt
end

# ╔═╡ dd57156e-3d76-470d-99d2-ded75692aaa4
reinforcement_models[last(findmax(M->last(M[2]),reinforcement_models))]

# ╔═╡ 9dba61da-6d76-4c9f-b3f1-fad8c2125b9b

let	(ship,path) = runModel(reinforcement_models[4,9][1], params)
	@info "remaining Δv" delta_v(last(path))
	@info "" apoapsis(last(path), params.body) periapsis(last(path), params.body)
	animate_path(
		path[begin:10:end], 
		params.body, 
		xlims=(0,20000),
		ylims=(0, 13100)
	)
end

# ╔═╡ dcf76654-bbe9-4234-b2c2-ce6c76cb08d7
gradient_model, gradient_rewards, gradient_path, gradient_grads = train_model!(
		θ -> ∇reward(θ, params),
		θ₀ = [0.5, 1.3, 1.22],
		steps=20,
		min_marginal_reward = 0.1m/s
	)

# ╔═╡ 34109d42-8bfe-4b99-8424-da28a9cf23bd
let	(ship,path) = runModel(gradient_model, params)
	@info "remaining Δv" delta_v(last(path))
	@info "" apoapsis(last(path), params.body) periapsis(last(path), params.body)
	animate_path(path[begin:end], params.body)
end

# ╔═╡ Cell order:
# ╠═6ddda634-6b63-11ed-0766-a57e14d8b741
# ╠═b051ae45-acc7-4de8-8ad6-abcff4fe6e17
# ╠═e2527ad1-1742-478e-bc52-7d9c20e535af
# ╠═102329cb-8c56-43f3-a1a4-7fab9984576a
# ╠═6c330af0-b2ea-4d8c-888f-b76356b38119
# ╠═1f0c5d25-cab5-4765-b0d5-473dc601778f
# ╠═c2290d22-e4eb-42c6-b01b-c13287caaab0
# ╠═08226305-1421-4a82-bb80-5df629f59957
# ╠═700255bf-88f8-4b72-9168-7bcfebf11184
# ╠═615bab31-67a4-4f56-923f-4011bafcba18
# ╠═ee3eb72f-741e-4e84-9809-346f30e9f84e
# ╠═8ba4a032-a56e-4d5c-8adf-fac41bf31d2c
# ╠═c88cbe4e-a4b6-4174-8d55-69efc4e8beaa
# ╠═30b8ca8e-d468-4afd-b0b7-f58b94b37246
# ╠═d11ab9f0-e8f9-46ee-85e4-596ff624cc65
# ╠═84dc7ae7-0f7c-45d4-9541-813f93348797
# ╠═04e5b976-58d5-4cfc-b95c-6bc52bbb1cfc
# ╠═89f65794-4041-40a1-a860-dedbef7fdf23
# ╠═1854c288-9503-4d99-82f6-b316db3070ef
# ╠═128a3180-5eee-4ffc-a53b-83833bb71ccc
# ╠═be0050c5-21eb-41ce-a2ea-18aaaf1d64a6
# ╠═c8be7b9c-4000-4f9f-a654-eace228b77cd
# ╠═58e0eb05-fe1b-4e59-bda6-b6f9bc95967d
# ╠═3f5a66b7-6b6d-4875-83be-8493e1d3724b
# ╠═03af5d6b-826a-45d9-8396-50a84e6bc07f
# ╠═856364cf-4804-4c5a-b2d3-b0c72697d90d
# ╠═d4f7bf69-ddb0-45d2-a895-6fb0794c5198
# ╠═bb53229e-2a85-4762-a1fb-980897d105d7
# ╠═71876131-4153-4fb9-97be-145922cadd4f
# ╠═de4659fc-2de4-44b6-bad9-4df78938f8d7
# ╠═0245f903-9feb-4b5f-abf3-3f34bb815dd4
# ╠═01905a71-be04-4214-a226-4d4a4686a1d2
# ╠═d5628d11-8301-4844-833d-3d77e4ca69dd
# ╠═0330467c-59a2-4941-b881-42fb0893e85b
# ╠═73ee0b41-8368-4f1f-9806-ef97780cbe7b
# ╠═2f621adb-7207-4831-8123-c081d66588a9
# ╠═3947e65c-ebd2-49ed-b842-4a2343457620
# ╠═537e3669-2861-4ac2-925c-a84a6f2d473a
# ╠═5dd533f3-b6c4-4591-ac3d-ec3960c6ece7
# ╠═efc466a6-5caf-4fdd-836f-c93c168b88ae
# ╠═c87cbaf2-4061-4991-9263-737509f2a945
# ╠═dd57156e-3d76-470d-99d2-ded75692aaa4
# ╠═9dba61da-6d76-4c9f-b3f1-fad8c2125b9b
# ╠═dcf76654-bbe9-4234-b2c2-ce6c76cb08d7
# ╠═34109d42-8bfe-4b99-8424-da28a9cf23bd
