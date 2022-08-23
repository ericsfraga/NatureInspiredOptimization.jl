# [[file:../book.org::tsppath][tsppath]]
struct Path
    city :: Vector{Int64}
    Path(n :: Int64) = new(1:n)
    Path(p :: Vector{Int64}) = new(copy(p))
end
Base.length(p :: Path) = length(p.city)
# tsppath ends here

# [[file:../book.org::tspindex][tspindex]]
Base.getindex(p :: Path, i :: Int64) = p.city[i]
function Base.setindex!(p :: Path, x, i :: Int64)
    p.city[i] = x
end
# tspindex ends here

# [[file:../book.org::tsp][tsp]]
function tsp(p :: Path, D :: Matrix{Float64})
    n = length(p)
    # sum of distances, closing the loop
    z = sum(D[p[i],p[i+1]] for i ∈ 1:n-1) + D[p[n],p[1]]
    # all paths are feasible
    g = 0.0
    (z, g)
end
# tsp ends here
