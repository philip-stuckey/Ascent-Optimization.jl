using FiniteDiff: finite_difference_gradient
using Unitful: m, s

∇reward(θ) = finite_difference_gradient(
	reward, 
	θ, 
	Val{:central},
	typeof(0.1m/s),
	relstep=0.001
)

function train_model!(∇reward; θ₀=zeros(3), steps=28, α=0.001s/m, min_marginal_reward=zero(RewardType))
	θ = θ₀[:]
	points = zeros(steps, length(θ))
	grads = zeros(RewardType, steps, length(θ))  # the fact that the gradient is technically the same type as the reward bothers me
	rewards = zeros(RewardType, steps)

	for n in 1:steps
		grad =  ∇reward(θ)
		θ .+= α .* grad
		r = reward(θ)
		if r - get(rewards,n-1, zero(r)) < min_marginal_reward
			θ .-= α .* grad
			@info "marginal reward too low" r-get(rewards,n-1,zero(r)) min_marginal_reward
            rewards = rewards[begin:n-1]
            points = points[begin:n-1, :]
            grads = grads[begin:n-1, :]
			break
		end
		rewards[n] = r
		points[n,:] .= θ
		grads[n,:] .= grad
	end
	return (model=Model(θ[1]/s,θ[2:end]), rewards=rewards, path=points, grads=grads)
end


