using Unitful
const 𝕃 = Unitful.𝐋
const 𝕋 = Unitful.𝐓

struct Planet
    gravitational_parameter::Quantity{Float64, 𝕃^3/𝕋^2}
    radius::Unitful.Length
    angular_speed::Unitful.Frequency
end

surface_gravity(body::Planet) = body.gravitational_parameter / body.radius^2

surface_speed(body::Planet) = body.radius * body.angular_speed

circular_orbit_speed(r, body) = √(body.gravitational_parameter/r)

