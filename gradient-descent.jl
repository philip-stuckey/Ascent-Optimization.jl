### A Pluto.jl notebook ###
# v0.19.14

using Markdown
using InteractiveUtils

# ╔═╡ d9e27bbc-5de6-11ed-3e0d-6d8758d91729
begin
	using Pkg
	Pkg.activate(".")
	include("src/OrbitalMechanics.jl")
end

# ╔═╡ e7b24203-e2a9-4b2f-8036-828a44a00085
begin
	using OrbitalMechanics
	using LinearAlgebra
end

# ╔═╡ f443c908-1ab1-4b2f-908e-43eb2bd7417b
begin
	using Plots
	gr()
	using Printf
end

# ╔═╡ cc1f1df3-a90f-4f7d-9859-be395eccf4e3
using FiniteDiff: finite_difference_gradient

# ╔═╡ 490b1fd9-2d88-406e-9dae-f7e1b629a0ff
body=Planet(1e8,1e3,0.0)

# ╔═╡ 58408b2d-d511-4744-8411-323b47919a2c
standard_parameters = OrbitalMechanics.SimulationParameters(
	body=body,
	initial_ship=Ship(body,Δv=2000, TWR=2, dry_mass=20),
	time_step=10^-5,
	margin=10.0,
	target_altitude=body.radius*1.5
)

# ╔═╡ 465c21ae-a3f3-4421-b8bf-9eb261ca6796
function loss(θ) 
	model = Model(θ[1], θ[2:end])
	rₘₐₓ = delta_v(standard_parameters.initial_ship)
	(ship, _) = runModel(model, standard_parameters, path=nothing)
	return (rₘₐₓ - delta_v(ship))/float(periapsis(ship,body) > standard_parameters.target_altitude)
end

# ╔═╡ a9fc7607-51bb-4d36-8559-b75f34d0f9a8
N=10

# ╔═╡ c1dbad33-c0af-4fcb-85f1-9a84b8f0371d
Ω = range(0.001,3, length=N)

# ╔═╡ 83b816c5-341a-4c93-82cd-a241990e0cee
θ₀ = range(-π/2, 0.3π, length=N)

# ╔═╡ 7cc4b431-6936-481b-9a3f-7e6927a3080a
θ₁ = range(-0.3, 0.44π, length=N)

# ╔═╡ 44f08774-9995-4dea-994d-8158b73bd641
# ╠═╡ disabled = true
#=╠═╡
coarse_search = [loss([ω, θ₀, θ₁]) for ω in Ω, θ₀ in θ₀, θ₁ in θ₁]
  ╠═╡ =#

# ╔═╡ f93b2938-7eae-4075-b132-98db3b8eecad
∇loss(v) = finite_difference_gradient(loss,v, relstep=0.1)

# ╔═╡ 4a2b752e-994a-40a6-96f4-472cfd7558d0
∇₀₀ = ∇loss(Float64[0.7, 0, 0])

# ╔═╡ d97aefcb-7e57-48ba-a493-f6924e2fec30
#=╠═╡
begin
	(u,v) = -normalize(∇₀₀) 
	
	heatmap(θ₀, θ₁, coarse_search[10,:,:], title="ω = $(Ω[10])", formatter=x-> @sprintf("%.2fπ", x/π))
	#scatter!([0],[0])
	#quiver!([0],[0], quiver=([u],[v]))
	#xlims!(extrema(θ₀)...)
	#ylims!(extrema(θ₁)...)
end
  ╠═╡ =#

# ╔═╡ b5020497-1014-4b62-a790-3a392afdeeac
function trainModel!(∇loss; ω₀=0.7, θ₀=zeros(3), steps=28, α=0.001)
	θ = vcat(ω₀, θ₀)
	points = zeros(steps,length(θ))
	grads = zeros(steps,length(θ))
	losses = zeros(steps)

	for n in 1:steps
		grad =  α .* ∇loss(θ)
		losses[n] = loss(θ)
		points[n,:] .= θ
		grads[n,:] .= grad
		θ .-= grad
	end
	return θ, points, losses, grads
end
	

# ╔═╡ b6da62f2-7d12-41ff-877c-9bef1a083ba0
trained_parameters, points, losses, grads = trainModel!(∇loss, steps=29, θ₀ = zeros(4))

# ╔═╡ ffeb4092-754b-4da3-b1c7-c4b1e8a2bc30
trained_model = Model(trained_parameters[1],trained_parameters[2:end])

# ╔═╡ bc1d469c-b18e-4a7c-90d6-e486f91a5725
begin 

	scatter(points[[1, end],1], points[[1,end],2], points[[1,end],3], color=[:cyan, :orange])
	plot!(points[:,1], points[:, 2], points[:, 3],legend=false)
	
	#quiver!(points[:,1], points[:,2], quiver=(grads[:,1], grads[:,2]))
	#lims!(extrema(θ₀)...)
	#ylims!(extrema(θ₁)...)
	#xlims!(0.6,0.75)
	#ylims!(-0.1,0.1)
end

# ╔═╡ 006e849b-3feb-4a99-8d3f-fb99eb2168fd
plot(losses, label="loss")

# ╔═╡ 86783c2f-a1ab-4390-b820-6e19010d8a55
2000 - last(losses)

# ╔═╡ d744f01f-2533-42aa-9d06-8c4edb361c10
let model = trained_model
	(ship,path) = runModel(model, standard_parameters)
	animate_path(path, body)
end

# ╔═╡ 74eb6000-c2cb-4afd-9ccc-a3a6d5779579
let model = trained_model
	(ship,path) = runModel(model, standard_parameters)
	plot([norm(ship.position) for ship in path[1:5]])
end

# ╔═╡ Cell order:
# ╠═d9e27bbc-5de6-11ed-3e0d-6d8758d91729
# ╠═e7b24203-e2a9-4b2f-8036-828a44a00085
# ╠═490b1fd9-2d88-406e-9dae-f7e1b629a0ff
# ╠═58408b2d-d511-4744-8411-323b47919a2c
# ╠═465c21ae-a3f3-4421-b8bf-9eb261ca6796
# ╠═a9fc7607-51bb-4d36-8559-b75f34d0f9a8
# ╠═c1dbad33-c0af-4fcb-85f1-9a84b8f0371d
# ╠═83b816c5-341a-4c93-82cd-a241990e0cee
# ╠═7cc4b431-6936-481b-9a3f-7e6927a3080a
# ╠═44f08774-9995-4dea-994d-8158b73bd641
# ╠═f443c908-1ab1-4b2f-908e-43eb2bd7417b
# ╠═d97aefcb-7e57-48ba-a493-f6924e2fec30
# ╠═cc1f1df3-a90f-4f7d-9859-be395eccf4e3
# ╠═f93b2938-7eae-4075-b132-98db3b8eecad
# ╠═4a2b752e-994a-40a6-96f4-472cfd7558d0
# ╠═b5020497-1014-4b62-a790-3a392afdeeac
# ╠═b6da62f2-7d12-41ff-877c-9bef1a083ba0
# ╠═ffeb4092-754b-4da3-b1c7-c4b1e8a2bc30
# ╠═bc1d469c-b18e-4a7c-90d6-e486f91a5725
# ╠═006e849b-3feb-4a99-8d3f-fb99eb2168fd
# ╠═86783c2f-a1ab-4390-b820-6e19010d8a55
# ╠═d744f01f-2533-42aa-9d06-8c4edb361c10
# ╠═74eb6000-c2cb-4afd-9ccc-a3a6d5779579