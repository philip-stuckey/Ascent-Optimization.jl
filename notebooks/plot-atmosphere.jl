### A Pluto.jl notebook ###
# v0.19.14

using Markdown
using InteractiveUtils

# ╔═╡ daa72a92-6512-11ed-0236-71db8e68e89b
begin
	import Pkg
	Pkg.activate("..")
	using AscentOptimization
end

# ╔═╡ ab4c5e99-319e-481b-90cb-007c3f0fd111
begin
	using Accessors
	using Unitful: kg,m,s
end

# ╔═╡ 0b9513be-c346-4ea8-9890-ab148b8c0f81
using LinearAlgebra, Plots

# ╔═╡ 511f941e-0b9c-4c5d-8df1-f569c944dd8f
parameters= let
	parameters= @set standard_parameters.body.atmosphere = AscentOptimization.ConstAtmosphere(1400m, 01.2kg/m^3)
	
end

# ╔═╡ 8b18523d-a188-42f8-8e40-2bb6afea7da9
plot_orbit(parameters.initial_ship, parameters.body)

# ╔═╡ 6a1450e4-9bd7-49fc-8d59-10d09b63c9e7
let model = Model(0.4/s, [0.0,0]), body = parameters.body
	(ship,path) = runModel(model, parameters)
	plot([specific_energy(ship,body) for ship in path])
end

# ╔═╡ d8c7c35f-38a3-4e35-98f0-292bed566b35
let model = Model(0.1/s, [0.0,0π]), body = parameters.body
	(ship,path) = runModel(model, parameters)
	animate_path(path, body)
end

# ╔═╡ b400850b-0884-473d-9101-000e87bcdbe5
Inf*oneunit(AscentOptimization.Length)

# ╔═╡ Cell order:
# ╠═daa72a92-6512-11ed-0236-71db8e68e89b
# ╠═ab4c5e99-319e-481b-90cb-007c3f0fd111
# ╠═511f941e-0b9c-4c5d-8df1-f569c944dd8f
# ╠═8b18523d-a188-42f8-8e40-2bb6afea7da9
# ╠═0b9513be-c346-4ea8-9890-ab148b8c0f81
# ╠═6a1450e4-9bd7-49fc-8d59-10d09b63c9e7
# ╠═d8c7c35f-38a3-4e35-98f0-292bed566b35
# ╠═b400850b-0884-473d-9101-000e87bcdbe5
