struct Planet
    gravitational_parameter::Float64
    radius::Float64
    angular_speed::Float64
end

g(body::Planet) = body.gravitational_parameter / body.radius^2