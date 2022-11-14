
using Plots
using Printf
using GeometryBasics: Point2
using Unitful: ustrip, unit

colTuple(v::Vector) = tuple(([u] for u in v)...)
function plot_orbit(ship, body; legend=false, aspect_ratio=1, kwargs...)
	
	ê = normalize(-eccentricity_vector(ship,body))
	ϕ = acos(ê[1:2]⋅[1,0])
	θ = range(0,2π, length=100)
	R = body.radius
	apoapsis = OrbitalMechanics.apoapsis(ship,body)
	periapsis = OrbitalMechanics.periapsis(ship,body)
	#o = orbit(ship, body)
	apoapsis_point = Point2(ê[1:2] * apoapsis)
	periapsis_point = Point2(-ê[1:2] *periapsis)
	
	plt = plot(;aspect_ratio, legend, kwargs...)
	# plot!(plt, -o.(θ) .* cos.(θ  .- ϕ), -o.(θ) .* sin.(θ .- ϕ))
	plot!(plt, R .* cos.(θ), R .* sin.(θ))
	scatter!(plt, [apoapsis_point, periapsis_point])
	scatter!(plt, Point2(ship.position[1:2]))
	annotate!(plt, [
		(apoapsis_point...," $(@sprintf("%.2f",ustrip(apoapsis)))$(unit(apoapsis)) "),
		(periapsis_point...," $(@sprintf("%.2f",ustrip(periapsis)))$(unit(periapsis)) ")
	])
	return plt
end

function animate_path(path, body)
	anim = @animate for (n,ship) in enumerate(path)
		plt1 = plot_orbit(ship,body)
		plot!(Point2.(ship.position[1:2] for ship in path[1:n]))
		quiver!(Point2(ship.position[1:2]), quiver=colTuple(100*heading(ship)[1:2]))
		e = eccentricity_vector(ship,body)*100
		quiver!([0],[0], quiver=([e[1], [e[2]]]))
		
		plt2 = plot(Point2.(ship.velocity[1:2] for ship in path[1:n]))
		scatter!(Point2(path[n].velocity[1:2]))

		plt3 = plot(Point2.(thrust_vector(ship)[1:2] for ship in path[1:n]))
		scatter!(Point2(thrust_vector(path[n])[1:2]))

		plot(
			plt1,
			#plt2, 
			#plt3
		)
	end 
	mov(anim)

end