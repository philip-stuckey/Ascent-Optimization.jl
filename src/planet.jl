struct Planet
    gravitational_parameter::Float64
    radius::Float64
    angular_speed::Float64
end

g(body::Planet) = body.gravitational_parameter / body.radius^2
circular_orbit_speed(r, body) = √(body.gravitational_parameter/r)
