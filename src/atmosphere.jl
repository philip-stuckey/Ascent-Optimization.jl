
abstract type AbstractAtmosphere end

using LinearAlgebra

function airspeed(ship,body)
    (_, τ) = basis(ship.position)
    v = ship.velocity
    r = norm(ship.position)
    ω = body.angular_speed
    return v - r*ω*τ
end

function drag_force(ship,body)
    height = norm(ship.position)
    ρ = density(body.atmosphere, height)
    v = airspeed(ship,body)
    D = drag_coefficient(ship)
    drag = ρ * (v⋅v) * D 
    return -v * drag
end

struct ConstAtmosphere <: AbstractAtmosphere
    height::Float64
    density::Float64
end

function density(atmosphere::ConstAtmosphere, height) 
    return atmosphere.density * (height > atmosphere.height)
end

using Interpolations
struct InterpolatedAtmosphere <: AbstractAtmosphere
    heights::Vector{Float64}
    density::Vector{Float64}
end

function density(atmosphere::InterpolatedAtmosphere, height)
    density_fnc = cubic_spline_interpolation(
        atmosphere.heights, 
        atmosphere.pressures
    )
    return density_fnc(height)
end
