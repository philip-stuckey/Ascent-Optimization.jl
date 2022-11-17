using Test
using .AscentOptimization: basis

η(ϕ) = Vec(sincos(ϕ)...,0)
@test all([ρ̂ ⋅ τ̂ == 0 for (ρ̂,τ̂) in basis.(η.(range(0,2π, length=50)))])