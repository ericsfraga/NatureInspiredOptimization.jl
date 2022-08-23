# [[file:../book.org::pfrmodelrhs][pfrmodelrhs]]
function reactor(y, u :: Real, x)
    dy = [ u * (15*y[2] - y[1])
           u * (y[1]-15*y[2])-(2-u)*y[2]]
    [ dy[1]
      dy[2]
      -dy[1]-dy[2]]
end
# pfrmodelrhs ends here

# [[file:../book.org::pfrsimulation][pfrsimulation]]
using DifferentialEquations
function simulate(u)
    y = [1.0, 0.0, 0.0] # initial conditions
    profile = [y]
    nx = length(u)
    δx = 2/nx
    for i ∈ 1:nx
        xspan = ((i-1)*δx, i*δx)
        prob = ODEProblem(reactor, y, xspan, u[i])
        ret = DifferentialEquations.solve(prob)
        y = ret.u[end]
        push!(profile, y)
    end
    # return last result and the full profile
    (y, profile) 
end
# pfrsimulation ends here

# [[file:../book.org::pfrobjectivefunction][pfrobjectivefunction]]
function f(u)
    # simulate the reactor with this rate
    solution, profile = simulate(u)
    # the objective is to maximise the production of P
    z = - solution[3]           # yield of y[3]
    (z, 0.0)                    # always feasible
end
# pfrobjectivefunction ends here

# [[file:../book.org::pfrcatalyststructure][pfrcatalyststructure]]
struct CatalystAmount
    α
    β
    γ
    x
    function CatalystAmount(r)
        new( [r[1], r[4]],      # α
             [r[2], r[5]],      # β
             [(1.0-r[1])*r[3],  # γ₁
              (1.0-r[4])*r[6]], # γ₂
             [r[7], 2.0-r[8]] ) # x
    end
end
# pfrcatalyststructure ends here

# [[file:../book.org::pfrcatalystbasis][pfrcatalystbasis]]
function a(x, β)
    (π/2 - atan(20*β*(2*x-1))) / π
end
function f₁(u, x)
    u.α[1] * a(x, u.β[1]) + u.γ[1]
end
function f₂(u, x)
    u.α[2] * a(x, u.β[2]) + u.γ[2]
end
# pfrcatalystbasis ends here

# [[file:../book.org::pfrcatalystindex][pfrcatalystindex]]
function Base.getindex(u :: CatalystAmount, x)
    if x ≤ u.x[1] && u.x[1] > eps()
        return f₁(u, x/u.x[1])
    elseif u.x[1] ≤ x ≤ u.x[2] && (u.x[2] - u.x[1]) > eps()
        return f₁(u, 1.0) + (f₂(u, 0.0) - f₁(u, 1.0))/(u.x[2]-u.x[1]) * (x - u.x[1])
    elseif u.x[2] ≤ x && u.x[2] < 2.0
        return f₂(u, (x-u.x[2])/(2.0-u.x[2]))
    else
        error("$x not in bounds [0,2]")
    end
end
# pfrcatalystindex ends here

# [[file:../book.org::pfrmodelaltrhs][pfrmodelaltrhs]]
function reactor(y, u :: CatalystAmount, x)
    dy = [ u[x] * (15*y[2] - y[1])
           u[x] * (y[1]-15*y[2])-(2-u[x])*y[2]]
    [ dy[1]
      dy[2]
      -dy[1]-dy[2]]
end
# pfrmodelaltrhs ends here

# [[file:../book.org::pfraltsimulation][pfraltsimulation]]
using DifferentialEquations
function simulate(u :: CatalystAmount)
    y = [1.0, 0.0, 0.0] # initial conditions
    profile = [y]
    intervals = [(0.0, u.x[1])
                 (u.x[1], u.x[2])
                 (u.x[2], 2.0)]
    for xspan ∈ intervals
        if xspan[2] - xspan[1] > eps()
            prob = ODEProblem(reactor, y, xspan, u)
            ret = DifferentialEquations.solve(prob) 
            y = ret.u[end]
            profile = vcat(profile, ret.u)
        end
    end
    (y, profile) # result and full profile
end
# pfraltsimulation ends here

# [[file:../book.org::pfraltobjectivefunction][pfraltobjectivefunction]]
function altf(u :: Vector{Float64})
    # simulate the reactor with this rate
    solution, profile = simulate(CatalystAmount(u))
    # the objective is to maximise the production of P
    z = - solution[3]           # yield of y[3]
    (z, 0.0)                    # always feasible
end
# pfraltobjectivefunction ends here
