module OrbitalMechanics

using StaticArrays
using LinearAlgebra

const Vec = SVector{3, Float64}
basis(r⃗) = (normalize(r⃗), [0 1 0; -1 0 0; 0 0 0]*normalize(r⃗) )

include("ship.jl")

end # module OrbitalMechanics
