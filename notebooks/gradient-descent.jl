### A Pluto.jl notebook ###
# v0.19.14

using Markdown
using InteractiveUtils

# ╔═╡ d9e27bbc-5de6-11ed-3e0d-6d8758d91729
begin
	using Pkg
	Pkg.activate("..")
end

# ╔═╡ 20d7a4cc-66f9-4345-a60f-f86d147edefe
begin
	using OrbitalMechanics
	using ThreadsX
	using Base.Iterators: product
end

# ╔═╡ ded14bf9-dcf1-432f-b21a-e0a8927a0c8d
using Unitful

# ╔═╡ 70954042-99b5-4937-9abf-b99997254594
using Base.Threads

# ╔═╡ bd927ee2-5f30-4865-83b4-0ef307c190b1
using Unitful: s,m

# ╔═╡ f443c908-1ab1-4b2f-908e-43eb2bd7417b
begin
	using Plots
	using Printf
	using GeometryBasics: Point2, Point3
end

# ╔═╡ 82cf3afa-c50a-47b4-b98e-2d5a05da3712
N=10

# ╔═╡ ab75e29b-2bdd-4641-b921-85402ed062fd
Ω = range(0.1,1.3, length=N)

# ╔═╡ 0833664f-d19c-40e5-bca3-ca8f526c0a14
Θ = range(0,π/2, length=N)

# ╔═╡ 7d7348ee-c73b-4f7c-bee0-b49bfe54b6b2
rewards = ThreadsX.map(x-> reward(Vec(x)), Iterators.product(Ω, Θ, Θ))

# ╔═╡ 801d53f1-e446-41f3-9b5d-a166d6c831c9
let M =1, θ₁ = Θ[M]
	heatmap(Ω, Θ, rewards[:,:,M], title="θ₂=$θ₁")
	#quiver!([0], [0], quiver=0.001 .* ([∇₀₀[2]], [∇₀₀[3]]))
end

# ╔═╡ b6da62f2-7d12-41ff-877c-9bef1a083ba0
let M = 3, steps=200, α=0.0001s/m, min_marginal_reward=0*m/s
	global models = []
	@threads for m in 3:M
		θ₀ = zeros(m)
		θ₀[1] = 0.2
		θ₀[end] = 2π
		push!(models, train_model!(∇reward; θ₀, steps, α, min_marginal_reward))
	end
end

# ╔═╡ 006e849b-3feb-4a99-8d3f-fb99eb2168fd
begin
	plot([m.rewards for m in models],legend=:bottomright)
	xlims!(0,110)
end

# ╔═╡ 8cf78cf9-1bbd-417f-a9b9-8b5f2ed441a4
plot(diff(models[1].rewards))

# ╔═╡ af736966-74ce-47aa-aa8f-d6386e87b4f9
let m = models[1].model
	reward(m)
end

# ╔═╡ 53fd6a63-21c9-400e-b99e-f61ab45aa7c9
let M =4, θ₁ = Θ[M]
	heatmap(Ω, Θ, rewards[:,:,M], title="θ₂=$θ₁")
	scatter!(Point2.(eachrow(models[1].path[:, [1,2]])))
end

# ╔═╡ d744f01f-2533-42aa-9d06-8c4edb361c10
let model = models[1].model
	(ship,path) = runModel(model, standard_parameters)
	animate_path(path, standard_parameters.body)
end

# ╔═╡ Cell order:
# ╠═d9e27bbc-5de6-11ed-3e0d-6d8758d91729
# ╠═82cf3afa-c50a-47b4-b98e-2d5a05da3712
# ╠═ab75e29b-2bdd-4641-b921-85402ed062fd
# ╠═0833664f-d19c-40e5-bca3-ca8f526c0a14
# ╠═20d7a4cc-66f9-4345-a60f-f86d147edefe
# ╠═7d7348ee-c73b-4f7c-bee0-b49bfe54b6b2
# ╠═801d53f1-e446-41f3-9b5d-a166d6c831c9
# ╠═70954042-99b5-4937-9abf-b99997254594
# ╠═b6da62f2-7d12-41ff-877c-9bef1a083ba0
# ╠═f443c908-1ab1-4b2f-908e-43eb2bd7417b
# ╠═006e849b-3feb-4a99-8d3f-fb99eb2168fd
# ╠═8cf78cf9-1bbd-417f-a9b9-8b5f2ed441a4
# ╠═af736966-74ce-47aa-aa8f-d6386e87b4f9
# ╠═53fd6a63-21c9-400e-b99e-f61ab45aa7c9
# ╠═d744f01f-2533-42aa-9d06-8c4edb361c10
