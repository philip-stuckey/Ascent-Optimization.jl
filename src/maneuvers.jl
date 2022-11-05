struct Maneuver
    done::Function
    throttle::Function
    declination::Function
end

function Maneuver(done, throttle::Float64=0, declination::Float64=0) 
    return Maneuver(done, (_...)-> throttle, (_...)-> declination)
end

function runManeuver!(ship, maneuver::Maneuver, parameters; path=nothing)
    body=parameters.body
    Δt=parameters.Δt
    snapshot_rate = parameters.snapshot_rate
	
    loops_since_last_snapshot=0
    for t in 0:Δt:parameters.time_limit
        maneuver.done(ship) && break
        delta_v(ship) <= 0 && break

        ship.throttle=maneuver.throttle(t)
        ship.declination = maneuver.declination(t)
        Simulate!(ship, body; Δt)
        loops_since_last_snapshot += 1

        if path !== nothing && loops_since_last_snapshot > 1/(snapshot_rate*Δt)
            push!(path, deepcopy(ship))
            loops_since_last_snapshot=0
        end
    end
    push!(path, deepcopy(ship))
    return path
end



end