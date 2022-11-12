using FiniteDiff: finite_difference_gradient

∇reward(θ) = finite_difference_gradient(reward, θ, relstep=0.01)

function trainModel!(∇reward; θ₀=zeros(3), steps=28, α=0.001, min_marginal_reward=0.0)
	θ = θ₀[:]
	points = zeros(steps,length(θ))
	grads = zeros(steps,length(θ))
	rewards = zeros(steps) .- 1

	for n in 1:steps
		grad =  ∇reward(θ)
		θ .+= α .* grad
		r = reward(θ)
		if r - get(rewards,n-1, 0) < min_marginal_reward
			θ .-= α .* grad
			@info "marginal reward too low" r-get(rewards,n-1,0) min_marginal_reward
			break
		end
		rewards[n] = r
		points[n,:] .= θ
		grads[n,:] .= grad
	end
	return (model=Model(θ[1],θ[2:end]), rewards=rewards, path=points, grads=grads)
end


