
abstract type AbstractAtmosphere end
using Unitful: kg, m
using LinearAlgebra

function airspeed(ship,body)
    (_, τ) = basis(ship.position)
    v = ship.velocity
    r = norm(ship.position)
    ω = body.angular_speed
    return v - r*ω*τ
end

function drag_force(ship, body)::Vec{Force}
    height = norm(ship.position)
    ρ = density(body.atmosphere, height)
    v = airspeed(ship,body)
    A = 2.0m^2  # drag coefficient
    d = 1.0
    drag = ρ * (v⋅v) * d * A 
    return (iszero(v) ? ustrip(v) : -normalize(v)) * drag
end

struct EmptyAtmosphere <: AbstractAtmosphere
end

density(::EmptyAtmosphere,r) = 0kg/m^3

struct ConstAtmosphere <: AbstractAtmosphere
    height::Length
    density::Density
end

function density(atmosphere::ConstAtmosphere, height) 
    return atmosphere.density * (height > atmosphere.height)
end

using Interpolations
struct InterpolatedAtmosphere <: AbstractAtmosphere
    heights::Vector{Length}
    density::Vector{Density}
end

function density(atmosphere::InterpolatedAtmosphere, height)
    density_fnc = cubic_spline_interpolation(
        atmosphere.heights, 
        atmosphere.pressures
    )
    return density_fnc(height)
end
