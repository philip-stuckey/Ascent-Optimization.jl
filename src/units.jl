using Unitful: Unitful, Quantity
using Unitful: m, s, kg
# using Unitful: Length, Mass, Time, Frequency, Velocity, Force

##=
const Length = Quantity{Float64, Unitful.𝐋, typeof(m)} 
const Mass = Quantity{Float64, Unitful.𝐌, typeof(kg)} 
const Time = Quantity{Float64, Unitful.𝐓, typeof(s)} 
const Frequency = Quantity{Float64, Unitful.𝐓^-1,typeof(s^-1)}

const Velocity = Quantity{Float64, Unitful.𝐋*Unitful.𝐓^-1, typeof(m*s^-1)} 
const Force = Quantity{Float64, Unitful.𝐌*Unitful.𝐋*Unitful.𝐓^-2, typeof(kg*m*s^-2)} 
# =#