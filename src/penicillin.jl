# [[file:../book.org::penicillinratestructure][penicillinratestructure]]
struct Rate
    U :: Vector{Float64}        # rate at each t point
    t :: Vector{Float64}        # t points
    function Rate(U, r, tf, δtmin)
        n = length(U)
        t = zeros(n)
        t[1] = 0.0
        for i ∈ 1:n-2
            δt = r[i]/2.0 * (tf - t[i])
            if δt < δtmin
                if t[i]+δtmin < tf
                    δt = δtmin
                end
            end
            t[i+1] = t[i] + δt
        end
        t[n] = tf
        new(U, t)
    end
end
# penicillinratestructure ends here

# [[file:../book.org::penicillinrateindex][penicillinrateindex]]
function Base.getindex(rate :: Rate, t :: Float64)
    n = length(rate.U)
    for i ∈ 2:n
        if t ≤ rate.t[i]
            return rate.U[i-1] +
                (t-rate.t[i-1]) *
                (rate.U[i] - rate.U[i-1]) /
                (rate.t[i] - rate.t[i-1])
        end
    end
    println("Warning: access time past end, $t > $(t[end])")
    rate.U[end]
end
Base.getindex(rate :: Rate, vt :: Vector{Float64}) = [rate[t] for t ∈ vt]
Base.firstindex(rate :: Rate) = rate.t[1]
Base.lastindex(rate :: Rate) = rate.t[end]
# penicillinrateindex ends here

# [[file:../book.org::penicillinmodelconstants][penicillinmodelconstants]]
μmax = 0.11                 # 1/h
ρmax = 0.0055               # g/g/h
Kx = 0.006                  # g/g
Kp = 0.0001                 # g/L
Kin = 0.1                   # g/L
Kdegr = 0.01                # 1/h
Km = 0.0001                 # g/L
ms = 0.029                  # g/g/h
Yxs = 0.47                  # g/g
Yps = 1.2                   # g/g
SF = 500                    # g/L
# penicillinmodelconstants ends here

# [[file:../book.org::penicillinreactormodel][penicillinreactormodel]]
function reactor(y, rate, t)
    # y = [X P S V]
    X, P, S, V = y
    U = rate[t]
    μ = μmax * S/(Kx*X + S)
    ρ = ρmax * S/(Kp + S*(1+S/Kin))
    [ μ*X-(X/(SF*V))*U
      ρ*X-Kdegr*P-(P/(SF*V))*U
      -μ*(X/Yxs)-ρ*(X/Yps)-(ms*S/(Km+S))*X+(1-S/SF)*U/V
      U/SF]
end
# penicillinreactormodel ends here

# [[file:../book.org::penicillinsimulation][penicillinsimulation]]
using DifferentialEquations
function simulate(rate :: Rate)
    y = [1.5  # X(0)
         0.0  # P(0)
         0.0  # S(0)
         7.0] # V(0)
    # also bounds on variables: X P S V
    upper = [40.0, 1e6, 100.0, 10.0]
    # used for constraint violations
    g = 0.0
    # save intermediate results for plotting
    profile = [y]
    lastt = 0.0
    for t ∈ rate.t
        # @show "loop start:", y, g
        if t>lastt
            tspan = (lastt, t)
            prob = ODEProblem(reactor, y, tspan, rate)
            ret = DifferentialEquations.solve(prob) # , maxiters=Int(1e6))
            y = ret.u[end]
            push!(profile, y)
            # @show "loop middle:", y, g
            lastt = t
            maxima = [maximum([ret.u[j][i] for j ∈ 1:length(ret.u)]) for i ∈ 1:4]
            violation = maximum(maxima - upper)
            if violation > g
                g = violation
            end
            # @show "loop end:", y, g
        end
    end
    (y, g, profile) # last result from solve and constraint
end
# penicillinsimulation ends here

# [[file:../book.org::penicillinobjectivefunction][penicillinobjectivefunction]]
function f₁(x :: Vector{Float64})
    # deconstruct argument:
    # x => U,   r, tf
    # len: n, n-2,  1
    #    = 2n-1
    n = Int((length(x)+1)/2)     # |U|
    U = x[1:n]
    r = x[n+1:end-1]
    tf = x[end]
    δtmin = tf/100
    # create rate structure
    rate = Rate(U, r, tf, δtmin)
    # simulate the reactor with this rate
    solution, g = simulate(rate)
    # the objective is to maximise the production of P
    z = - solution[2] * solution[4]
    (z, g)
end
# penicillinobjectivefunction ends here
