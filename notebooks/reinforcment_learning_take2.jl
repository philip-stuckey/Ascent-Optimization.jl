### A Pluto.jl notebook ###
# v0.19.14

using Markdown
using InteractiveUtils

# ╔═╡ f63f90fd-4666-416e-b294-365582e17ae8
begin
	using Pkg;
	Pkg.activate("..")
	using Interpolations
	using Plots
	using LinearAlgebra
	using AscentOptimization
	using Printf
	using GeometryBasics: Point2
	using Base.Threads
	using ThreadsX
	md"# Deriving Optimal Ascent Trajectory with Reinforcement Learning"
end

# ╔═╡ da64679a-1374-467d-a27c-490661f6903b
using Unitful: kg, m, s

# ╔═╡ 67b5f6ff-8916-4bfa-afcc-ab787de6c063
let model=Model(declination=[0,π/2+0.08], pitch_rate=0.4/s)
	N = 1:5
	rewards = AscentOptimization.Velocity[]
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
			time_step=(10.0^-n)s
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

# ╔═╡ 6d8e3c80-be6f-4f14-9712-3785e3f9be18
εs=range(0.1,1.0, length=10)

# ╔═╡ ea8f64e8-205b-43e2-a78c-bac7229620c8
αs=range(0.1, 3, length=10)

# ╔═╡ 69e7e51b-0930-43bf-8d9b-f2901cca780a
begin
	global models = fill(Model(0.0/s, []), length(εs), length(αs))
	global rewards = zeros(AscentOptimization.Velocity, length(εs), length(αs))
	for (n,α) in collect(enumerate(αs))
		for (m,ε) in enumerate(εs)
			k=1
			# for k in 1:100
				step_rewards = []
				model =train_model!(
					Model(declination=[0, π/2], pitch_rate=1.0/s),
					standard_parameters,
					EpsilonExplorer(α, ε),
					steps=100, 
					rewards=step_rewards
				)
				models[n,m] = model
				rewards[n,m] += (last(step_rewards) - rewards[n,m])/k
			# end
		end
	end
end

# ╔═╡ 94be906b-2de3-4b58-afa3-e2bf79fb703f
heatmap(αs, εs, rewards)

# ╔═╡ bd793b4b-9990-429c-be6a-47318a3073d1
trained_model, training_rewards = let explorer = EpsilonExplorer(0.4,0.5)
	model = Model(declination=[0, π/2], pitch_rate=1.0/s)
	step_rewards = AscentOptimization.RewardType[]
	model =train_model!(model, standard_parameters, explorer, steps=100, rewards=step_rewards)
	
	(model, step_rewards)
end

# ╔═╡ 87491a89-4b87-4fc9-844a-360a0736db6e
random_model = let	
	global random_model_rewards = AscentOptimization.RewardType[]
	model =train_model!(
		Model(declination=[0, π/2], pitch_rate=1.0/s), 
		standard_parameters, 
		EpsilonExplorer(0.4,1),
		steps=1000, 
		rewards=random_model_rewards
	)
end

# ╔═╡ d1ac734b-a637-43b9-ba61-12433bb0b24c
trained_model.declination ./ π

# ╔═╡ 38079515-96db-4454-af3b-63830da3b6d0
last(training_rewards)

# ╔═╡ 5c26ae3f-0e93-43dd-bb5f-4c27bbd9ee7c
last(random_model_rewards)

# ╔═╡ ecc57a7c-9e67-4d17-95c1-76bc1ad62dd3
begin
	plot(legend=:bottomright)
	plot!(random_model_rewards, label="random explorer")
	plot!(training_rewards, label="ε=0.8")
end

# ╔═╡ 5695e0cf-1e51-4624-af95-0d1e3c34314f
let (_, path) = runModel(trained_model, standard_parameters)
	@info reward(trained_model, standard_parameters)
	animate_path(path, standard_parameters.body)
end

# ╔═╡ 5d4abc9a-0e40-495d-8ce9-074cdd43caa2
let (_, path) = runModel(random_model, standard_parameters)
	@info reward(random_model, standard_parameters)
	animate_path(path, standard_parameters.body)
end

# ╔═╡ Cell order:
# ╠═f63f90fd-4666-416e-b294-365582e17ae8
# ╠═da64679a-1374-467d-a27c-490661f6903b
# ╠═67b5f6ff-8916-4bfa-afcc-ab787de6c063
# ╠═6d8e3c80-be6f-4f14-9712-3785e3f9be18
# ╠═ea8f64e8-205b-43e2-a78c-bac7229620c8
# ╠═69e7e51b-0930-43bf-8d9b-f2901cca780a
# ╠═94be906b-2de3-4b58-afa3-e2bf79fb703f
# ╠═bd793b4b-9990-429c-be6a-47318a3073d1
# ╠═87491a89-4b87-4fc9-844a-360a0736db6e
# ╠═d1ac734b-a637-43b9-ba61-12433bb0b24c
# ╠═38079515-96db-4454-af3b-63830da3b6d0
# ╠═5c26ae3f-0e93-43dd-bb5f-4c27bbd9ee7c
# ╠═ecc57a7c-9e67-4d17-95c1-76bc1ad62dd3
# ╠═5695e0cf-1e51-4624-af95-0d1e3c34314f
# ╠═5d4abc9a-0e40-495d-8ce9-074cdd43caa2
