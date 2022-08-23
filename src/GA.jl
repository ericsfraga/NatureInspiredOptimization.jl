# [[file:../book.org::crossover][crossover]]
function crossover(p1 :: Vector{T}, p2 :: Vector{T}) where T <: Number
    n = length(p1)
    @assert n == length(p2) "Points to crossover must be of same length"
    # generate vector of true/false values with
    # true meaning "selected" in comments that
    # follow
    alleles = rand(Bool, n)
    # first new solution (n1) consists of p1 with
    # selected alleles from p2; note that we need
    # to copy the arguments as otherwise we would
    # be modifying the contents of the arguments.
    n1 = copy(p1)
    n1[alleles] = p2[alleles]
    # the second new solution is the other way
    # around: p2 with selected from p1
    n2 = copy(p2)
    n2[alleles] = p1[alleles]
    n1, n2
end
# crossover ends here

# [[file:../book.org::mutate][mutate]]
function mutate(p, lower, upper)
    m = copy(p)      # passed by reference
    n = length(p)    # vector length
    i = rand(1:n)    # index for value to mutate
    a = lower[i]     # lower bound for that value
    b = upper[i]     # and upper bound
    r = rand()       # r ∈ [0,1]
    # want r ∈ [-0.5,0.5] to mutate around p[i]
    m[i] = p[i] + (b-a) * (r-0.5)
    # check for boundary violations and rein back in
    m[i] = m[i] < a ? a : (m[i] > b ? b : m[i])
    m
end
# mutate ends here

# [[file:../book.org::ga][ga]]
using NatureInspiredOptimization: fitness, mostfit, Point, select, statusoutput
function ga(
    # required arguments, in order
    p0,                   # initial population
    f;                    # objective function
    # optional arguments, in any order
    lower = nothing,      # bounds for search domain
    upper = nothing,      # bounds for search domain
    parameters = nothing, # for objective function 
    nfmax = 10000,        # max number of evaluations
    np = 100,             # population size
    output = true,        # output during evolution
    rc = 0.7,             # crossover rate
    rm = 0.05)            # mutation rate

    p = p0                # starting population
    lastmag = 0           # for status output 
    lastnf = 0            # also for status output
    nf = 0                # count evaluations
    while nf < nfmax
        ϕ = fitness(p)
        best = mostfit(p,ϕ)
        newp = [best]       # elite size 1
        lastmag, lastnf = statusoutput(output, nf,
                                       best, lastmag,
                                       lastnf)
        i = 0
        while i < np
            print(stderr, "nf=$nf i=$i\r")
            point1 = p[select(ϕ)]
            if rand() < rc      # crossover?
                point2 = p[select(ϕ)]
                newpoint1, newpoint2 = crossover(point1.x, point2.x)
                push!(newp,
                      Point(newpoint1, f, parameters))
                push!(newp,
                      Point(newpoint2, f, parameters))
                nf += 2
                i += 2
            elseif rand() < rm  # mutate?
                newpoint1 = mutate(point1.x, lower, upper)
                push!(newp,
                      Point(newpoint1, f, parameters))
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
    lastmag, lastnf = statusoutput(output, nf,
                                   best,
                                   lastmag, lastnf)
    best, p, ϕ
end
# ga ends here
