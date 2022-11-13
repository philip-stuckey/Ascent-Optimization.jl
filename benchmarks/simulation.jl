using OrbitalMechanics
using BenchmarkTools
using Unitful: s

ship = standard_parameters.initial_ship
body = standard_parameters.body
Δt=0.001s

@info "Simulate!" @benchmark OrbitalMechanics.Simulate!($(deepcopy(ship)), $body; Δt = $Δt)

@info "Eulars!" @benchmark OrbitalMechanics.Eulars!($(deepcopy(ship)), $body; Δt = $Δt) 

@info "Stick to Ground" @benchmark OrbitalMechanics.stick_to_ground!($(deepcopy(ship)), $body, $Δt) 