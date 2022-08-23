# [[file:../book.org::pso][pso]]
using NatureInspiredOptimization: ≻, Point, randompoint, statusoutput
function pso(
    # required arguments, in order
    p0,                   # initial population
    f;                    # objective function
    # optional arguments, in any order
    lower = nothing,      # bounds for search domain
    upper = nothing,      # bounds for search domain
    parameters = nothing, # for objective function 
    nfmax = 10000,        # max number of evaluations
    np = 10,              # population size
    output = true,        # output during evolution
    c₁ = 2.0,             # weight for velocity
    c₂ = 2.0,             # ditto
    ω = 1.0               # inertia for velocity
    )

    # create initial full population including p0
    p = copy(p0)
    for i ∈ length(p0)+1:np
        push!(p, Point(randompoint(lower, upper),
                       f, parameters))
    end
    nf = np
    # find current best in population
    best = begin
        b = p[1]
        for i ∈ 2:np
            if p[i] ≻ b
                b = p[i]
            end
        end
        b
    end
    h = p            # historical best
    v = []           # initial random velocity
    for i ∈ 1:np
        push!(v, (upper .- lower) .*
            ( 2 * rand(length(best.x)) .- 1))
    end
    lastmag = 0      # for status output
    lastnf = 0       # as well
    while nf < nfmax
        lastmag, lastnf = statusoutput(output, nf,
                                       best,
                                       lastmag,
                                       lastnf)
        for i ∈ 1:np
            print(stderr, "nf=$nf i=$i\r")
            v[i] = ω * v[i] .+
                c₁ * rand() * (h[i].x .- p[i].x) .+
                c₂ * rand() * (best.x .- p[i].x)
            x = p[i].x .+ v[i]
            # fix any bounds violations
            x[x .< lower] = lower[x .< lower]
            x[x .> upper] = upper[x .> upper]
            # evaluate and save new position
            p[i] = Point(x, f, parameters)
            nf += 1
            # update history for this point
            if p[i] ≻ h[i]
                h[i] = p[i]
            end
            # and global best found so far
            if p[i] ≻ best
                best = p[i]
            end
        end
    end
    lastmag, lastnf = statusoutput(output, nf, best,
                                   lastmag, lastnf)
    best, p
end
# pso ends here
