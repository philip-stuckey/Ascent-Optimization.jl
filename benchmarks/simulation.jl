using AscentOptimization
using BenchmarkTools
using Unitful: s

ship = standard_parameters.initial_ship
body = standard_parameters.body
Δt=0.001s

@info "Simulate!" @benchmark AscentOptimization.Simulate!($(deepcopy(ship)), $body; Δt = $Δt)

@info "Eulars!" @benchmark AscentOptimization.Eulars!($(deepcopy(ship)), $body; Δt = $Δt) 

@info "Stick to Ground" @benchmark AscentOptimization.stick_to_ground!($(deepcopy(ship)), $body, $Δt) 