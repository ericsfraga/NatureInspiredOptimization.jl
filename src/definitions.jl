# [[file:../book.org::succfloat][succfloat]]
≻(a :: Number, b :: Number) = a < b
≽(a :: Number, b :: Number) = a ≤ b
export ≻, ≽
# succfloat ends here

# [[file:../book.org::dominates][dominates]]
dominates(a, b) = all(a .≽ b) && any(a .≻ b)
# dominates ends here

# [[file:../book.org::succvector][succvector]]
≻(a :: Vector{T}, b :: Vector{T}) where T <: Number = dominates(a,b)
# succvector ends here

# [[file:../book.org::mostfit][mostfit]]
mostfit(p,ϕ) = p[ϕ .≥ maximum(ϕ)][1]
# mostfit ends here

# [[file:../book.org::nondominated][nondominated]]
function nondominated(pop)
    nondom = Point[]
    for p1 ∈ pop
        dominated = false
        for p2 ∈ pop
            if p2 ≻ p1
                dominated = true
                break
            end
        end
        if ! dominated
            push!(nondom, p1)
        end
    end
    nondom
end
# nondominated ends here

# [[file:../book.org::printpoints][printpoints]]
function printpoints(name, points)
    println("#+name: $name")
    for p ∈ points
        printvector("| ", p.x, "| ")
        printvector("| ", p.z, "| ")
        printvector("| ", p.g, "| ")
        println()
    end
end
# printpoints ends here

# [[file:../book.org::printvector][printvector]]
function printvector(start, x, separator = "")
    @printf "%s" start
    for i ∈ 1:length(x)
        @printf "%.3f %s" x[i] separator
    end
end
# printvector ends here

# [[file:../book.org::randompoint][randompoint]]
randompoint(a, b) = a + (b - a) .* rand(length(a))
# randompoint ends here

# [[file:../book.org::statusoutput][statusoutput]]
function statusoutput(
    output,  # true/false
    nf,      # number of function evaluations
    best,    # population of Points
    lastmag, # magnitude of last output
    lastnf   # nf at last output
    )
    if output
        z = best.z
        if length(z) == 1
            z = z[1]
        end
        δ = nf - lastnf
        mag = floor(log10(nf))
        if mag > lastmag
            println("evolution: $nf $z $(best.g)")
            lastmag = mag
            lastnf = nf
        elseif δ > 10^(mag)-1
            println("evolution: $nf $z $(best.g)")
            lastnf = nf
        end
    end
    (lastmag, lastnf)
end
# statusoutput ends here

# [[file:../book.org::point][point]]
struct Point
    x :: Any       # decision point
    z :: Vector    # objective function values
    g :: Real      # constraint violation
    function Point(x, f, parameters = nothing)
        # evaluate the objective function,
        # using parameters if given
        if parameters isa Nothing
            z, g = f(x)
        else
            z, g = f(x, parameters)
        end
        # check types of returned values
        # and convert if necessary
        if g isa Int
            g = float(g)
        end
        if rank(z) == 1
            # vector of objective function values
            new(x, z, g)
        elseif rank(z) == 0
            # if scalar, create size 1 vector
            new(x, [z], g)
        else
            error( "NatureInspiredOptimization" *
                " methods can only handle scalar" *
                " and vector criteria," *
                " not $(typeof(z))." )
        end
    end
end
# point ends here

# [[file:../book.org::rank][rank]]
rank(x :: Any) = length(size(x))
# rank ends here

# [[file:../book.org::succpoint][succpoint]]
≻(a :: Point, b :: Point) = (a.g ≤ 0) ? (b.g > 0.0 || a.z ≻ b.z) : (a.g < b.g)
# succpoint ends here

# [[file:../book.org::fitness][fitness]]
function fitness(pop :: Vector{Point})
    l = length(pop)
    # feasible solutions are those with g ≤ 0
    indexfeasible = (1:l)[map(p->p.g,pop) .<= 0]
    # infeasible have g > 0
    indexinfeasible = (1:l)[map(p->p.g,pop) .> 0]
    fit = zeros(l)
    # the factor will be used to squeeze the
    # fitness values returned by vectorfitness into
    # the top half of the interval (0,1) for
    # feasible solutions and the bottom half for
    # infeasible, if both types of solutions are
    # present.  Otherwise, the full interval is
    # used.
    factor = 1   
    if length(indexfeasible) > 0
        # consider only feasible subset of pop
        feasible = view(pop,indexfeasible)
        # use objective function value(s) for ranking
        feasiblefit = vectorfitness(map(p->p.z,feasible))
        if length(indexinfeasible) > 0
            # upper half of fitness interval
            feasiblefit = feasiblefit./2 .+ 0.5
            # have both feasible & infeasible
            factor = 2 
        end
        fit[indexfeasible] = (feasiblefit.+factor.-1)./factor
    end
    if length(indexinfeasible) > 0
        # squeeze infeasible fitness values into
        # (0,0.5) or (0,1) depending on factor,
        # i.e. whether there are any feasible
        # solutions as well or not
        infeasible = view(pop,indexinfeasible)
        # use constraint violation for ranking as
        # objective function values although it
        # should be noted that the measure of
        # infeasibility may not actually make any
        # sense given that points are infeasible
        fit[indexinfeasible] = vectorfitness(map(p->p.g, infeasible)) / factor;
    end
    fit
end
# fitness ends here

# [[file:../book.org::vectorfitness][vectorfitness]]
function vectorfitness(v)
    # determine number of objectives (or
    # pseudo-objectives) to consider in ranking
    l = length(v)
    if l == 1
        # no point in doing much as there is only
        # one solution
        [0.5]
    else
        m = length(v[1])      # number of objectives
        if m == 1             # single objective 
            fitness = [v[i][1] for i=1:l]
        else                  # multi-objective
            # rank of each solution for each
            # objective function
            rank = ones(m,l);
            for i=1:m
                rank[i,sortperm([v[j][i] for j=1:l])] = 1:l
            end
            # hadamard product of ranks
            fitness = map(x->prod(x), rank[:,i] for i=1:l)
        end
        # normalise (1=best, 0=worst) while
        # avoiding extreme 0,1 values using the
        # hyperbolic tangent
        adjustfitness(fitness)
    end
end
# vectorfitness ends here

# [[file:../book.org::adjustfitness][adjustfitness]]
function adjustfitness(fitness)
    if (maximum(fitness)-minimum(fitness)) > eps()
        0.5*(tanh.(4*(maximum(fitness) .- fitness) /
            (maximum(fitness)-minimum(fitness)) .- 2)
             .+ 1)
    else
        # if there is only one solution in the
        # population or if all the solutions are
        # the same, the fitness value is halfway
        # in the range allowed
        0.5 * ones(length(fitness))
    end
end
# adjustfitness ends here

# [[file:../book.org::select][select]]
function select(ϕ)
    n = length(ϕ)
    @assert n > 0 "Population is empty"
    if n ≤ 1
        1
    else
        i = rand(1:n)
        j = rand(1:n)
        ϕ[i] > ϕ[j] ? i : j
    end
end
# select ends here

# [[file:../book.org::hypervolume][hypervolume]]
function hypervolume(points)
    m = length(points)
    n = length(points[1].z)
    set = Matrix{Real}(undef, m, n)
    for r ∈ 1:m
        set[r,1:n] = points[r].z
    end
    z = sortslices(set, dims=1, by=x->(x[1], x[2]))
    area = 0.0
    last = z[m,n]
    first = z[1,1]
    distance = sqrt(z[1,1]^2 + z[1,2]^2)
    for i ∈ 2:m
        area += (z[i,1]-z[i-1,1]) * ((z[1,2]-z[i-1,2])+(z[1,2]-z[i,2])) / 2.0
        d = sqrt(z[i,1]^2 + z[i,2]^2)
        if d < distance
            distance = d
        end
    end
    # maximise area, minimise distance
    # larger outcome better
    area/distance^2
end
# hypervolume ends here

# [[file:../book.org::randomsearch][randomsearch]]
using NatureInspiredOptimization: Point, statusoutput
function randomsearch(
    # required arguments, in order
    f;                    # objective function
    # optional arguments, in any order
    lower = nothing,      # bounds for search domain
    upper = nothing,      # bounds for search domain
    parameters = nothing, # for objective function 
    nfmax = 10000,        # max number of evaluations
    output = true)        # output during evolution

    lastmag = 0           # for status output 
    lastnf = 0            # also for status output

    best = Point(randompoint(lower,upper), f,
                 parameters)
    nf = 1                # function evaluations
    while nf < nfmax
        p = Point(randompoint(lower,upper),
                  f, parameters)
        nf += 1
        if p ≻ best
            best = p
        end
        lastmag, lastnf = statusoutput(output, nf,
                                       best,
                                       lastmag,
                                       lastnf)
    end
    lastmag, lastnf = statusoutput(output, nf,
                                   best,
                                   lastmag, lastnf)
    best 
end
# randomsearch ends here
