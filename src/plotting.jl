
using Plots
using Printf

colTuple(v::Vector) = tuple(([u] for u in v)...)
function plot_orbit(ship, body; legend=false, aspect_ratio=1, kwargs...)
	
	ê = normalize(eccentricity_vector(ship,body))
	ϕ = acos(ê[1:2]⋅[1,0])
	θ = range(0,2π, length=100)
	R = body.radius
	#o = orbit(ship, body)
	apoapsis_point = Point2(ê[1:2] * apoapsis(ship,body))
	periapsis_point = Point2(-ê[1:2] * periapsis(ship,body))
	
	plt = plot(;aspect_ratio, legend, kwargs...)
	# plot!(plt, -o.(θ) .* cos.(θ  .- ϕ), -o.(θ) .* sin.(θ .- ϕ))
	plot!(plt, R .* cos.(θ), R .* sin.(θ))
	
	scatter!(plt, [apoapsis_point, periapsis_point])
	scatter!(plt, Point2(ship.position[1:2]))
	annotate!(plt, [
		(apoapsis_point...,"$(@sprintf("%.2f",apoapsis(ship,body)))"),
		(periapsis_point...,"$(@sprintf("%.2f",periapsis(ship,body)))")
	])
	return plt
end