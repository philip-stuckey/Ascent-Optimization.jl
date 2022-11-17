### A Pluto.jl notebook ###
# v0.19.14

using Markdown
using InteractiveUtils

# ╔═╡ ca1df62e-66ac-11ed-0818-c9be4223555c
begin
	using Pkg
	Pkg.activate("..")
end

# ╔═╡ 64fd12d9-bf17-4dbf-9203-487359f704a5
begin
	using AscentOptimization
	using ThreadsX
	using Base.Iterators: product
end

# ╔═╡ 56219bd3-cb59-41b8-99ce-12b1365d5477
using Plots: heatmap

# ╔═╡ 269533c9-5823-4ebf-bd8b-f76bb179a0ad
N=10

# ╔═╡ 554b6480-6451-4b10-87d1-b4a27d64e69b
Ω = range(0.1,1.3, length=N)

# ╔═╡ 146d0e84-fd91-4e40-ac8e-69108d141258
Θ = range(0,π/2, length=N)

# ╔═╡ 8f80e040-add2-4ea1-82b7-f916b4cde9ad
rewards = ThreadsX.map(x-> reward(Vec(x)), Iterators.product(Ω, Θ, Θ));

# ╔═╡ 89d13997-4e24-440f-9985-764d70dc5265
let M =3
	heatmap(Ω, Θ, rewards[M,:,:], title="ω=$(Ω[M])", xlabel="θ₁", ylabel="θ₂")
end

# ╔═╡ 5b3928c0-8d05-4e48-9e4b-b2ce20bdb8d2
let M =3
	heatmap(Ω, Θ, rewards[:,M,:], title="θ₁=$(Θ[M])", xlabel="ω", ylabel="θ₂")
end

# ╔═╡ f3070aa4-bb65-4cd8-a957-45f815418cca
let M =3
	heatmap(Ω, Θ, rewards[:,:,M], title="θ₂=$(Θ[M])", xlabel="ω", ylabel="θ₁")
end

# ╔═╡ Cell order:
# ╠═ca1df62e-66ac-11ed-0818-c9be4223555c
# ╠═269533c9-5823-4ebf-bd8b-f76bb179a0ad
# ╠═554b6480-6451-4b10-87d1-b4a27d64e69b
# ╠═146d0e84-fd91-4e40-ac8e-69108d141258
# ╠═64fd12d9-bf17-4dbf-9203-487359f704a5
# ╠═8f80e040-add2-4ea1-82b7-f916b4cde9ad
# ╠═56219bd3-cb59-41b8-99ce-12b1365d5477
# ╠═89d13997-4e24-440f-9985-764d70dc5265
# ╠═5b3928c0-8d05-4e48-9e4b-b2ce20bdb8d2
# ╠═f3070aa4-bb65-4cd8-a957-45f815418cca
