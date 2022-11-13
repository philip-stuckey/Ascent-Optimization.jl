using Unitful
using Unitful: m, s, kg
using Unitful: Length, Mass, Time, Frequency, Velocity, Force

#=
const Length = Quantity{Float64, Unitful.𝐋, m}
const Mass = Quantity{Float64, Unitful.𝐌, kg}
const Time = Quantity{Float64, Unitful.𝐓, s}
const Frequency = Quantity{Float64, Unitful.𝐓^-1, s^-1}
const Velocity = Quantity{Float64, Unitful.𝐋*Unitful.𝐓^-1, m*s^-1}
const Force = Quantity{Float64, Unitful.𝐌*Unitful.𝐋*Unitful.𝐓^-2, kg*m*s^-2}
=#