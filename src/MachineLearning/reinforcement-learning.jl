using LinearAlgebra
using Base.Iterators: product, flatten

randargmax(A) = rand([n for (n,a) in enumerate(A) if a .== maximum(A)])

struct EpsilonExplorer
    learning_rate::Float64
    exploration_rate::Float64
end

struct QModel{A}
	quality::Dict{A, Float64}
    count::Dict{A, Int}
end

QModel{R}(actions) where R = QModel(Dict(actions.=>zero(R)), Dict(actions.=>0))

function choice(explorer::EpsilonExplorer, model::QModel)
	actions = collect(keys(model.quality))
	vals = values(model.quality)
	return if explorer.learning_rate >= rand()
		rand(actions)
	else
		actions[randargmax(vals)]
	end
end

function update!(model::QModel, action, reward, learning_rate)
	model.count[action] += 1
	model.quality[action] += learning_rate * (reward - model.quality[action])
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

function train_model!(
	m₀::Model, 
	simulation_parameters::SimulationParameters, 
	explorer::EpsilonExplorer; 
	rewards=nothing,
	steps=100
)
	r₀ = reward(m₀, simulation_parameters)
	
	learner = QModel(action_space(m₀))
	model = m₀
	for _ in 1:steps
		rewards == nothing || push!(rewards, r₀)
		action = choice(explorer, learner)
		model = apply(model, action)
		r = reward(model, simulation_parameters)
		update!(learner, action, r-r₀, explorer.learning_rate)
		if r < r₀
			model = apply(model, negate_action(action))
		else
			r₀ = r
		end
	end
	return model
end