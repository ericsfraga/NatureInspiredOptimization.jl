# [[file:../book.org::package][package]]
module NatureInspiredOptimization

# external packages required
using Printf # for the @printf macro

# utility functions and common structures for the nature inspired
# optimization methods
include("definitions.jl")

# chlorobenzene process optimization
module Chlorobenzene
include("chlorobenzene.jl")
end

# penicillin fermentation process operation
module PFR
include("PFR.jl")
end

# heat exchanger network design
module HEN
include("HEN.jl")
end

# the genetic algorithm
module GA
include("GA.jl")
end

# the particle swarm optimization method
module PSO
include("PSO.jl")
end

# the plant propagation algorithm
module PPA
include("PPA.jl")
end

end # module
# package ends here
