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

# ╔═╡ 70954042-99b5-4937-9abf-b99997254594
using Base.Threads

# ╔═╡ f443c908-1ab1-4b2f-908e-43eb2bd7417b
begin
	using Plots
	using Printf
	using GeometryBasics: Point2, Point3
end

# ╔═╡ ded14bf9-dcf1-432f-b21a-e0a8927a0c8d
using Unitful

# ╔═╡ bd927ee2-5f30-4865-83b4-0ef307c190b1
using Unitful: s,m

# ╔═╡ 82cf3afa-c50a-47b4-b98e-2d5a05da3712
N=10

# ╔═╡ ab75e29b-2bdd-4641-b921-85402ed062fd
Ω = range(0.1,1.3, length=N)

# ╔═╡ 0833664f-d19c-40e5-bca3-ca8f526c0a14
Θ = range(0,π/2, length=N)

# ╔═╡ 7d7348ee-c73b-4f7c-bee0-b49bfe54b6b2
rewards = ThreadsX.map(x-> reward(Vec(x)), Iterators.product(Ω, Θ, Θ))

# ╔═╡ e3bd314c-1744-49bf-a1b7-6e3c7eb663a5
reward(Vec((0.1, 0, 0)))

# ╔═╡ 7122065f-02b3-414d-81af-b7411b537c8c
∇reward([0.1,0,0])

# ╔═╡ c23c4b2f-e095-4bf9-934d-34eefa5541cd
ω, θ₁, θ₂ = (0.1, 0, π/2)

# ╔═╡ 3d3a5674-5989-4bb0-8723-03ec72cbefc2
δω, δθ₁, δθ₂ = ∇reward([ω, θ₁, θ₂]) .|> ustrip

# ╔═╡ 801d53f1-e446-41f3-9b5d-a166d6c831c9
let M =1, θ = Θ[M]
	heatmap(Ω, Θ, rewards[:,:, M], title="θ₂=$θ")
	quiver!([ω], [θ], quiver=([δω/1000], [δθ₁/1000]))
	scatter!([ω], [θ])
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
	#xlims!(0,110)
end

# ╔═╡ 8cf78cf9-1bbd-417f-a9b9-8b5f2ed441a4
plot(diff(models[1].rewards))

# ╔═╡ af736966-74ce-47aa-aa8f-d6386e87b4f9
let m = models[1].model
	reward(m)
end

# ╔═╡ 53fd6a63-21c9-400e-b99e-f61ab45aa7c9
let M =4, θ₁ = Θ[M], N = 
	heatmap(Ω, Θ, rewards[:,:,M], title="final Δv where θ₂=$(@sprintf("%.2f",θ₁))")

	(model, _, path, grads) = models[1]
	
	ωs, θ₁s, θ₂s = eachcol(path)
	scatter!(ωs[[1,end]], θ₁s[[1,end]])
	plot!(ωs, θ₁s)
	xlabel!("pitch rate")
	ylabel!("initial pitch")
	
	#δωs, δθ₁s, δθ₂s = eachcol(grads[[1,end],:]) ./ 100000 .|> ustrip
	#quiver!(ωs, θ₁s, quiver=(δωs, δθ₁s))
	#xlims!(0.0, 0.5)
end

# ╔═╡ d744f01f-2533-42aa-9d06-8c4edb361c10
let model = models[1].model, body = standard_parameters.body
	(ship,path) = runModel(model, standard_parameters)
	animate_path(path,body)
end

# ╔═╡ Cell order:
# ╠═d9e27bbc-5de6-11ed-3e0d-6d8758d91729
# ╠═82cf3afa-c50a-47b4-b98e-2d5a05da3712
# ╠═ab75e29b-2bdd-4641-b921-85402ed062fd
# ╠═0833664f-d19c-40e5-bca3-ca8f526c0a14
# ╠═20d7a4cc-66f9-4345-a60f-f86d147edefe
# ╠═7d7348ee-c73b-4f7c-bee0-b49bfe54b6b2
# ╠═e3bd314c-1744-49bf-a1b7-6e3c7eb663a5
# ╠═7122065f-02b3-414d-81af-b7411b537c8c
# ╠═c23c4b2f-e095-4bf9-934d-34eefa5541cd
# ╠═ded14bf9-dcf1-432f-b21a-e0a8927a0c8d
# ╠═3d3a5674-5989-4bb0-8723-03ec72cbefc2
# ╠═801d53f1-e446-41f3-9b5d-a166d6c831c9
# ╠═70954042-99b5-4937-9abf-b99997254594
# ╠═bd927ee2-5f30-4865-83b4-0ef307c190b1
# ╠═b6da62f2-7d12-41ff-877c-9bef1a083ba0
# ╠═f443c908-1ab1-4b2f-908e-43eb2bd7417b
# ╠═006e849b-3feb-4a99-8d3f-fb99eb2168fd
# ╠═8cf78cf9-1bbd-417f-a9b9-8b5f2ed441a4
# ╠═af736966-74ce-47aa-aa8f-d6386e87b4f9
# ╠═53fd6a63-21c9-400e-b99e-f61ab45aa7c9
# ╠═d744f01f-2533-42aa-9d06-8c4edb361c10
# ╠═ded14bf9-dcf1-432f-b21a-e0a8927a0c8d
# ╠═bd927ee2-5f30-4865-83b4-0ef307c190b1
