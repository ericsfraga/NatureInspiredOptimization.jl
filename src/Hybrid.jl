# [[file:../book.org::steepestdescent][steepestdescent]]
using NatureInspiredOptimization: printvector
function steepestdescent(
    # required arguments
    x0,                         # initial point
    f;                          # objective function
    # optional arguments
    lower = nothing,            # lower bounds
    upper = nothing,            # upper bounds
    parameters = nothing,       # for objective
    ϵ = 1e-4                    # for convergence test
    )
    # keep track of number of function evaluations
    nf = 0
    # function for line search:
    # α: step
    # x: starting point
    # g: gradient at x
    # c: 1 for function value, 2 for constraint violation
    h(α,x,g) = (parameters isa Nothing ? f(x - α * g) : f(x - α * g, parameters))
    x = x0
    while true                  # iterate
        # check feasibility of current point
        z = parameters isa Nothing ? f(x) : f(x,parameters)
        c = z[2] ≤ 0.0 ? 1 : 2     # feasible?
        g = ∇(f, x, parameters, c) # find grad
        # perform the line search starting with x
        ( α, nnf ) = search(x, g, h, c, lower, upper, ϵ)
        nf += nnf                           # keep count
        # printvector("x: ", x)
        # print(" f: $(f(x))")
        # printvector(" g: ", g)
        # println(" α = $α")
        if α > ϵ
            x = x - α * g       # update point
        else
            # return result and function evaluations taken
            return (x, nf)
        end
    end                         # exit guaranteed
end
# steepestdescent ends here

# [[file:../book.org::grad][grad]]
# unit vector in i-th dimension
unitvector(i,n) = [j == i for j in 1:n]
# f: objective function
# x: point
# p: parameters, if necessary
# c: objective function value or constraint violation
function ∇(f, x, p, c)
    δ = 1e-8
    n = length(x)
    # if infeasible, no checks required
    if c == 2
        if p isa Nothing
            [(f(x.+δ*unitvector(i,n))[2] - f(x)[2])/δ for i in 1:n]
        else
            [(f(x.+δ*unitvector(i,n),p)[2] - f(x,p)[2])/δ for i in 1:n]
        end 
    else
        g = []
        z = p isa Nothing ? f(x) : f(x,p)
        for i in 1:n
            if p isa Nothing
                z2 = f(x .+ δ * unitvector(i,n))
            else
                z2 = f(x .+ δ * unitvector(i,n), p)
            end
            if z2[2] > 0
                # infeasible; check left difference
                if p isa Nothing
                    z2 = f(x .- δ * unitvector(i,n))
                else
                    z2 = f(x .- δ * unitvector(i,n), p)
                end
                if z2[2] > 0
                    # must lie on equality constraint: do not move
                    # from here
                    push!(g, 0.0)
                else
                    push!(g, (z[1]-z2[1])/δ)
                end
            else
                push!(g, (z2[1]-z[1])/δ)
            end
        end
        g
    end
end
# grad ends here

# [[file:../book.org::search][search]]
function search(x, g, h, c, lower, upper, ϵ)
    α = 1.0
    h0 = h(0.0, x, g)
    z0 = h0[c]
    # println("search: c=$c x=$x h=$h0")
    nf = 0
    while α > ϵ
        # check that bounds are not violated
        if ( lower isa Nothing || all(x.+α*g .≥ lower) ) &&
            ( upper isa Nothing || all(x.+α*g .≤ upper))
            # take the step
            h2 = h(α, x, g)
            nf += 1
            # println("α=$α h2=$h2")
            # if feasible, new point must also be feasible
            # println("test: $(((c == 1 && h2[2] ≤ 0)))")
            # println("test: $(((c == 1 && h2[2] ≤ 0) || c == 2))")
            # println("test: $(((c == 1 && h2[2] ≤ 0) || c == 2) && h2[c] < z0)")
            if ((c == 1 && h2[2] ≤ 0) || c == 2) && h2[c] < z0
                # found a better point
                # println("search: h($α)=$(h(α,x,g))")
                return ( α, nf )
            end
        end
        # bisect
        α = α/2.0
    end
    # better solution not found
    return ( 0.0, nf )
end
# search ends here

# [[file:../book.org::add!][add!]]
function add!(pop, point, ϵ=1e-6)
    if all([maximum(abs.(p.x-point.x))>ϵ for p in pop])
        push!(pop, point)
    end
end
# add! ends here

# [[file:../book.org::hybridga][hybridga]]
using NatureInspiredOptimization: fitness, mostfit, Point, select, statusoutput
using NatureInspiredOptimization.GA: crossover, mutate
using NatureInspiredOptimization.Hybrid: add!, steepestdescent
function hybridga(
    # required arguments, in order
    p0,                     # initial population
    f;                      # objective function
    # optional arguments, in any order
    lower = nothing,        # lower bounds for search domain
    upper = nothing,        # upper bounds for search domain
    parameters = nothing,   # extra parameters for objective function 
    nfmax = 10000,          # maximum number of function evaluations
    np = 100,               # population size
    output = true,          # print output during evolution
    rc = 0.7,               # crossover rate
    rm = 0.05,              # mutation rate
    ϵ = 0.001               # tolerance for diversity control
    )        

    p = p0                  # starting population
    lastmag = 0             # for status output during evolution
    lastnf = 0              # also for status output
    nf = 0                  # count of function evaluations
    while nf < nfmax
        ϕ = fitness(p)      # fitness of current population
        best = mostfit(p,ϕ)
        newp = [best]       # elite size 1
        # apply local procedure to best
        ( localpoint, nnf ) = steepestdescent(best.x, f;
                                              lower = lower,
                                              upper = upper,
                                              parameters = parameters,
                                              ϵ = ϵ)
        nf += nnf           # keep count of function evaluations
        add!(newp, Point(localpoint, f, parameters), ϵ)
        lastmag, lastnf = statusoutput(output, nf, best,
                                       lastmag, lastnf)
        i = 0
        while i < np
            print(stderr, "nf=$nf i=$i\r")
            point1 = p[select(ϕ)]
            if rand() < rc      # crossover?
                point2 = p[select(ϕ)]
                newpoint1, newpoint2 = crossover(point1.x, point2.x)
                push!(newp, Point(newpoint1, f, parameters))
                push!(newp, Point(newpoint2, f, parameters))
                nf += 2
                i += 2
            elseif rand() < rm  # mutate?
                newpoint1 = mutate(point1.x, lower, upper)
                push!(newp, Point(newpoint1, f, parameters))
                nf += 1
                i += 1
            else                # copy over
                push!(newp, point1) 
                i += 1
            end
        end
        p = newp
    end
    ϕ = fitness(p)
    best = mostfit(p,ϕ)
    lastmag, lastnf = statusoutput(output, nf, best, lastmag, lastnf)
    best, p, ϕ
end
# hybridga ends here

# [[file:../book.org::hybridpso][hybridpso]]
using NatureInspiredOptimization: ≻, Point, randompoint, statusoutput
function hybridpso(
    # required arguments, in order
    p0,                         # initial population
    f;                          # objective function
    # optional arguments, in any order
    lower = nothing,        # lower bounds for search domain
    upper = nothing,        # upper bounds for search domain
    parameters = nothing,   # extra parameters for objective function 
    nfmax = 10_000,         # maximum number of function evaluations
    np = 10,                # population size
    output = true,          # print output during evolution
    c₁ = 2.0,               # as per original paper
    c₂ = 2.0,               # ditto
    ϵ = 1e-3,               # tolerance for local solver
    ω = 0.7                 # inertia for velocity change
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
        push!(v, (upper .- lower) .* ( 2 * rand(length(best.x)) .- 1))
    end
    lastmag = 0      # for status output during evolution
    lastnf = 0       # also for status output
    while nf < nfmax
        lastmag, lastnf = statusoutput(output, nf, best,
                                       lastmag, lastnf)
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
                # if we have found better point,
                # apply local procedure and then
                # save
                ( localpoint, nnf ) = steepestdescent( p[i].x, f;
                                                       lower = lower,
                                                       upper = upper,
                                                       parameters = parameters,
                                                       ϵ = ϵ)
                nf += nnf           # keep count of function evaluations
                best = Point(localpoint, f, parameters)
            end
        end
    end
    lastmag, lastnf = statusoutput(output, nf, best,
                                   lastmag, lastnf)
    best, p
end
# hybridpso ends here

# [[file:../book.org::hybridppa][hybridppa]]
using NatureInspiredOptimization: fitness, mostfit, Point, select, statusoutput
using NatureInspiredOptimization.PPA: neighbour
using NatureInspiredOptimization.Hybrid: add!, steepestdescent
function hybridppa(
    # required arguments, in order
    p0,                     # initial population
    f;                      # objective function
    # optional arguments, in any order
    lower = nothing,        # lower bounds for search domain
    nfmax = 10000,          # maximum number of function evaluations
    np = 10,                # population size
    nrmax = 5,              # maximum number of runners
    output = true,          # print output during evolution
    parameters = nothing,   # extra parameters for objective function 
    ϵ = 0.001,              # tolerance for diversity control
    upper = nothing         # upper bounds for search domain
    )

    p = p0                  # starting population
    lastmag = 0             # for status output during evolution
    lastnf = 0              # also for status output
    nf = 0                  # count of function evaluations
    while nf < nfmax
        ϕ = fitness(p)      # fitness of current population
        best = mostfit(p,ϕ)
        lastmag, lastnf = statusoutput(output, nf, best, lastmag, lastnf)
        newp = [best]       # elite set size 1
        # apply local search procedure
        # @localsolve steepestdescent f best parameters lower upper newp 1e-6 1e-6
        ( localpoint, nnf ) = steepestdescent(best.x, f;
                                              lower = lower,
                                              upper = upper,
                                              parameters = parameters,
                                              ϵ = ϵ)
        nf += nnf       # keep count of function evaluations
        add!(newp, Point(localpoint, f, parameters), ϵ)
        for i ∈ 1:np
            print(stderr, "nf=$nf i=$i\r")
            s = select(ϕ)
            nr = ceil(rand() * ϕ[s] * nrmax)
            for j ∈ 1:nr
                newpoint = neighbour(p[s].x, lower, upper, ϕ[s])
                add!(newp, Point(newpoint, f, parameters), ϵ)
                nf += 1
            end
        end
        p = newp
    end
    ϕ = fitness(p)
    best = mostfit(p,ϕ)
    lastmag, lastnf = statusoutput(output, nf, best, lastmag, lastnf)
    best, p, ϕ
end
# hybridppa ends here
