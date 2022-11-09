### A Pluto.jl notebook ###
# v0.19.14

using Markdown
using InteractiveUtils

# ╔═╡ f63f90fd-4666-416e-b294-365582e17ae8
begin
	using Pkg;
	Pkg.activate(".")
	using Interpolations
	using Plots
	using LinearAlgebra
	include("src/OrbitalMechanics.jl")
	using OrbitalMechanics
	using Printf
	using GeometryBasics: Point2
	using Base.Threads
	using ThreadsX
	md"# Deriving Optimal Ascent Trajectory with Reinforcement Learning"
end

# ╔═╡ 648ac1b0-cf45-46cc-959c-e7f3bbb5a51f
begin
	struct EpsilonExplorer{A}
		learning_rate::Float64
		exploration_rate::Float64
		quality::Dict{A, Float64}
		N::Dict{A, Int}
	end
	EpsilonExplorer(ε, α, actions) = EpsilonExplorer(ε, α, Dict(actions .=> 0.0), Dict(actions .=> 0))
end

# ╔═╡ 8c0b77b5-e538-4549-bb41-c23b9364b686
randargmax(A) = rand([n for (n,a) in enumerate(A) if a .== maximum(A)])

# ╔═╡ d302fbc3-9315-45e2-9e9b-c711a0516b32
function choice(explorer::EpsilonExplorer)
	actions = collect(keys(explorer.quality))
	vals = values(explorer.quality)
	return if explorer.learning_rate >= rand()
		rand(actions)
	else
		actions[randargmax(vals)]
	end
end

# ╔═╡ a95486a8-ec44-4ef9-9f07-1b9e5dcc5b39
function update!(explorer::EpsilonExplorer, action, reward)
	explorer.N[action] += 1
	explorer.quality[action] += explorer.learning_rate * (reward - explorer.quality[action])
end

# ╔═╡ 7297bae3-afa9-46d4-aa1f-062eece787a1
function action_space(model::Model)
	n = length(model.declination)
	basis = eachcol(Matrix(I, n,n))
	Δθ = deg2rad(1)
	
	pitch_space = Iterators.flatten((-1,1) .* Ref(Δθ * ê) for ê in basis)
	rate_space = [-0.1, 0.1]
	return collect(Iterators.product(pitch_space, rate_space))[:]
end

# ╔═╡ dc4c5959-79eb-4475-b212-f3e589d522c8
action_space(Model(declination=[0,π/2+0.08], pitch_rate=0.4)) 

# ╔═╡ a624f6f2-02a5-4dcb-89ac-958ad3b5efa6
function apply(model::Model, action)
	(Δθ, Δω) = action
	θ₀ = model.declination
	ω₀ = model.pitch_rate
	return Model(declination = θ₀ .+ Δθ, pitch_rate = ω₀ + Δω)
end

# ╔═╡ 7125b4c6-84eb-43a3-b2b2-522ac7efcdf7
negate_action(action) = -1 .* action

# ╔═╡ f0fa75ff-45e2-47c1-bc6c-36b2f653c0b2
Base.@kwdef struct TrainingParameters
	learning_rate::Float64
	exploration_rate::Float64
	step_limit::Int
	max_reward::Float64=Inf
end

# ╔═╡ 43bf9daf-531e-405f-b93b-34042e46dc75
function reward(model::Model, params::SimulationParameters)
	(ship, _) = runModel(model, params, path=nothing)
	return delta_v(ship)
end

# ╔═╡ 67b5f6ff-8916-4bfa-afcc-ab787de6c063
let model=Model(declination=[0,π/2+0.08], pitch_rate=0.4)
	N = 1:5
	rewards = []
	times = []

	# warm up
	reward(model, standard_parameters) 
	for n in N
		params = SimulationParameters(
			body=standard_parameters.body, 
			initial_ship=standard_parameters.initial_ship, 
			target_altitude=standard_parameters.target_altitude,
			margin=standard_parameters.margin,
			time_limit=standard_parameters.time_limit,
			time_step=10.0^-n
		)
		res = @timed reward(model, params) 
		push!(rewards, res.value)
		push!(times, res.time)
	end
	reward_plot = plot(N,rewards, xlabel="-log(Δt)", ylabel="reward")
	time_plot = plot(N, times, xlabel="-log(Δt)", ylabel="simulation time (seconds)")
	cross_plot = plot(times, rewards, xlabel="simulation time (seconds)", ylabel="reward")
	plot(reward_plot, time_plot, cross_plot, legend=false)
end

# ╔═╡ dcd28735-8ce4-4717-8b1a-11ca6b1c2d80
function train_model!(
	m₀::Model, 
	simulation_parameters::SimulationParameters, 
	params::TrainingParameters; 
	rewards=nothing
)
	r₀ = reward(m₀, simulation_parameters)
	explorer = EpsilonExplorer(
		params.learning_rate, 
		params.exploration_rate, 
		action_space(m₀)
	)
	
	model = m₀
	for n in 1:params.step_limit
		rewards == nothing || push!(rewards, r₀)
		action = choice(explorer)
		model = apply(model, action)
		r = reward(model, simulation_parameters)
		update!(explorer, action, r-r₀)
		if r < r₀
			model = apply(model, negate_action(action))
		else
			r₀ = r
		end
		#@info "got here" rewards r rewards==nothing params.max_reward
		
		r > params.max_reward && break
	end
	return model
end

# ╔═╡ bd793b4b-9990-429c-be6a-47318a3073d1
trained_model, training_rewards = let params = TrainingParameters(
	learning_rate=0.8,
	exploration_rate=0.5,
	step_limit=200,
	max_reward=Inf
)
	model = Model(declination=[0, π/2], pitch_rate=1.0)
	step_rewards = []
	model =train_model!(model, standard_parameters, params, rewards=step_rewards)
	
	(model, step_rewards)
end

# ╔═╡ d1ac734b-a637-43b9-ba61-12433bb0b24c
trained_model.declination ./ π

# ╔═╡ 38079515-96db-4454-af3b-63830da3b6d0
last(training_rewards)

# ╔═╡ 8f672f7d-39f5-4ece-931b-aaaac1b0c36b
plot(training_rewards,legend=false)

# ╔═╡ 5695e0cf-1e51-4624-af95-0d1e3c34314f
let (_, path) = runModel(trained_model, standard_parameters)
	animate_path(path, standard_parameters.body)
end

# ╔═╡ ecc57a7c-9e67-4d17-95c1-76bc1ad62dd3
let	params = TrainingParameters(
	learning_rate=0.4,
	exploration_rate=1,
	step_limit=200,
	max_reward=Inf
)
	model = Model(declination=[0, π/2], pitch_rate=1.0)
	step_rewards = []
	model =train_model!(model, standard_parameters, params, rewards=step_rewards)
	
	plot(step_rewards, label=last(step_rewards), title="random explorer")
end

# ╔═╡ Cell order:
# ╠═f63f90fd-4666-416e-b294-365582e17ae8
# ╟─67b5f6ff-8916-4bfa-afcc-ab787de6c063
# ╠═648ac1b0-cf45-46cc-959c-e7f3bbb5a51f
# ╠═8c0b77b5-e538-4549-bb41-c23b9364b686
# ╠═d302fbc3-9315-45e2-9e9b-c711a0516b32
# ╠═a95486a8-ec44-4ef9-9f07-1b9e5dcc5b39
# ╠═7297bae3-afa9-46d4-aa1f-062eece787a1
# ╠═dc4c5959-79eb-4475-b212-f3e589d522c8
# ╠═a624f6f2-02a5-4dcb-89ac-958ad3b5efa6
# ╠═7125b4c6-84eb-43a3-b2b2-522ac7efcdf7
# ╠═f0fa75ff-45e2-47c1-bc6c-36b2f653c0b2
# ╠═43bf9daf-531e-405f-b93b-34042e46dc75
# ╠═dcd28735-8ce4-4717-8b1a-11ca6b1c2d80
# ╠═bd793b4b-9990-429c-be6a-47318a3073d1
# ╠═d1ac734b-a637-43b9-ba61-12433bb0b24c
# ╠═38079515-96db-4454-af3b-63830da3b6d0
# ╠═8f672f7d-39f5-4ece-931b-aaaac1b0c36b
# ╠═5695e0cf-1e51-4624-af95-0d1e3c34314f
# ╠═ecc57a7c-9e67-4d17-95c1-76bc1ad62dd3
