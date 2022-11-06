using FiniteDiff: finite_difference_gradient

function loss(θ) 
    model = Model(θ[1], θ[2:end])
    rₘₐₓ = delta_v(standard_parameters.initial_ship)
    (ship, _) = runModel(model, standard_parameters, path=nothing)
    return (rₘₐₓ - delta_v(ship))/float(periapsis(ship,body) > standard_parameters.target_altitude)
end


