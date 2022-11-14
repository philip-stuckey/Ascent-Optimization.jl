using LinearAlgebra
using Base.Iterators: product, flatten

randargmax(A) = rand([n for (n,a) in enumerate(A) if a .== maximum(A)])

struct EpsilonExplorer{A,R}
    learning_rate::Float64
    exploration_rate::Float64
    quality::Dict{A, R}
    count::Dict{A, Int}
end

function EpsilonExplorer{R}(ε, α, actions) where R <: Number
    default_quality = Dict(actions .=> zero(R))
    default_count = Dict(actions .=> 0)
    return EpsilonExplorer(ε, α, default_quality, default_count)
end

function choice(explorer::EpsilonExplorer)
	actions = collect(keys(explorer.quality))
	vals = values(explorer.quality)
	return if explorer.learning_rate >= rand()
		rand(actions)
	else
		actions[randargmax(vals)]
	end
end

function update!(explorer::EpsilonExplorer, action, reward)
	explorer.N[action] += 1
	explorer.quality[action] += explorer.learning_rate * (reward - explorer.quality[action])
end


function action_space(model::Model)
	n = length(model.declination)
	basis = eachcol(Matrix(I, n,n))
	Δθ = deg2rad(1)
	
	pitch_space = flatten((-1,1) .* Ref(Δθ * ê) for ê in basis)
	rate_space = [-0.1, 0.1]
	return collect(product(pitch_space, rate_space))[:]
end

negate_action(action) = -1 .* action

function apply(model::Model, action)
	(Δθ, Δω) = action
	θ₀ = model.declination
	ω₀ = model.pitch_rate
	return Model(declination = θ₀ .+ Δθ, pitch_rate = ω₀ + Δω)
end

Base.@kwdef struct TrainingParameters
	learning_rate::Float64
	exploration_rate::Float64
	step_limit::Int
	max_reward::Float64
end

function train_model!(
	m₀::Model, 
	simulation_parameters::SimulationParameters, 
	params::TrainingParameters; 
	rewards=nothing
)
	r₀ = reward(m₀, simulation_parameters)
	explorer = EpsilonExplorer{Float64}(
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
        		
		r > params.max_reward && break
	end
	return model
end