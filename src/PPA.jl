# [[file:../book.org::neighbour][neighbour]]
function neighbour(s,     # selected point
                   lower, # lower bounds
                   upper, # upper bounds
                   ϕ)     # fitness
    r = rand(length(s))
    # element by element calculation
    ν = s + (upper - lower) .*
        ((2 * (r .- 0.5)) * (1-ϕ))
    # check for bounds violations
    ν[ν .< lower] = lower[ν .< lower]
    ν[ν .> upper] = upper[ν .> upper]
    return ν
end
# neighbour ends here

# [[file:../book.org::ppa][ppa]]
using NatureInspiredOptimization: fitness, mostfit, Point, select, statusoutput
function ppa(
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
    nrmax = 5)            # maximum number of runners

    p = p0                # starting population
    lastmag = 0           # for status output 
    lastnf = 0            # also for status output
    nf = 0                # function evaluations
    while nf < nfmax
        ϕ = fitness(p)
        best = mostfit(p,ϕ)
        newp = [best]       # elite set size 1
        lastmag, lastnf = statusoutput(output, nf,
                                       best,
                                       lastmag,
                                       lastnf)
        newpoints = [] 
        for i ∈ 1:np
            print(stderr, "nf=$nf i=$i\r")
            s = select(ϕ)
            nr = ceil(rand() * ϕ[s] * nrmax)
            for j ∈ 1:nr
                push!(newpoints,
                      neighbour(p[s].x, lower,
                                upper, ϕ[s]))
                nf += 1
            end
        end
        Threads.@threads for x in newpoints
            push!(newp, Point(x, f, parameters))
        end 
        p = newp
    end
    ϕ = fitness(p)
    best = mostfit(p,ϕ)
    lastmag, lastnf = statusoutput(output, nf,
                                   best,
                                   lastmag, lastnf)
    best, p, ϕ
end
# ppa ends here
