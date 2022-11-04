struct Planet
    gravitational_parameter::Float64
    radius::Float64
    angular_speed::Float64
end

surface_gravity(body::Planet) = body.gravitational_parameter / body.radius^2

surface_speed(body::Planet) = body.radius * body.angular_speed

circular_orbit_speed(r, body) = √(body.gravitational_parameter/r)

