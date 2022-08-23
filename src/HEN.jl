# [[file:../book.org::henstreamstructure][henstreamstructure]]
mutable struct Stream
    name   # descriptive
    type   # :hot or :cold
    mcp    # W or kW / K
    Tin    # K or C
    Tout   # K or C
    h      # W or kW / K / m^2
    Q      # W or kW typically
end
# henstreamstructure ends here

# [[file:../book.org::*Stream][Stream:2]]
# constructor which calculates Q
Stream(name,type,mcp,Tin,Tout,h) =
    Stream(name, type, mcp,
           Tin, Tout, h, 
           abs(mcp*(Tout-Tin)))
# constructor which copies existing Stream instance
Stream(s :: Stream) = Stream(s.name, s.type, s.mcp,
                             s.Tin, s.Tout, s.h, s.Q)
# Stream:2 ends here

# [[file:../book.org::*Stream][Stream:3]]
# informative output
Base.show(io :: IO, s :: Stream) = print(io, "$(s.name): $(s.Tin) → $(s.Q)$(s.type == :hot ? "↓" : "↑") → $(s.Tout)")
# predicates
iscold(s :: Stream) = s.type == :cold
ishot(s :: Stream) = s.type == :hot
# Stream:3 ends here

# [[file:../book.org::henutilitystructure][henutilitystructure]]
struct ExternalUtility
    name
    type           # hot or cold
    Tin
    Tout
    h
    model          # f(Q)
end
ishot(u :: ExternalUtility) = u.type == :hot
# henutilitystructure ends here

# [[file:../book.org::heninfeasibility][heninfeasibility]]
struct Infeasibility <: Exception
    g  # measure of infeasibility
end
# heninfeasibility ends here

# [[file:../book.org::henexchanger][henexchanger]]
struct Exchanger
    hot     # hot stream temperatures
    cold    # cold stream temperatures
    ΔT      # log-mean temperature difference
    U       # overall heat transfer coefficient
    A       # area [m²]
    Q       # duty [kW]
    cost    # capital cost in currency units
end
# henexchanger ends here

# [[file:../book.org::*Heat exchanger <<henheatexchanger>>][Heat exchanger <<henheatexchanger>>:2]]
function Exchanger(hot, cold, Q, model)
    # define tuple of temperatures: (in, out)
    if hot isa Stream
        h = (hot.Tin, hot.Tin-Q/hot.mcp)
        if h[2] < hot.Tout
            throw(Infeasibility(hot.Tout-h[2]))
        end
    else
        h = (hot.Tin, hot.Tout)
    end
    # (in, out) for cold side as well
    if cold isa Stream
        c = (cold.Tin, cold.Tin+Q/cold.mcp)
        if c[2] > cold.Tout
            throw(Infeasibility(c[2]-cold.Tout))
        end
    else
        c = (cold.Tin, cold.Tout)
    end
    # temperature driving forces
    ΔT1 = h[1]-c[2]
    ΔT2 = h[2]-c[1]
    if ΔT1 ≤ 0.0 || ΔT2 ≤ 0.0 || Q ≤ 0.0
        # indicate infeasibility using duty as
        # measure
        throw(Infeasibility(abs(Q)))
    end
    # Chen approximation to log mean temperature
    # difference
    ΔT = ( ΔT1 * ΔT2 * (ΔT1 + ΔT2)/2.0) ^ (1.0/3.0)
    # overall heat transfer coefficient
    U = 1.0 / (1.0/hot.h + 1.0/cold.h)
    # Area
    A = Q/U/ΔT
    # capital cost as a function of the area
    cost = model(A) # 6600.0 + 670.0 * A^0.83
    # create the object
    Exchanger(h, c, ΔT, U, A, Q, cost)
end
# Heat exchanger <<henheatexchanger>>:2 ends here

# [[file:../book.org::*Heat exchanger <<henheatexchanger>>][Heat exchanger <<henheatexchanger>>:3]]
using Printf # make @printf available
function Base.show(io :: IO, e :: Exchanger)
    @printf(io,
            "HEX(H:%.1f-%.1f C:%.1f-%.1f Q=%.2f: ΔT=%.2f U=%.3f A=%.2f %.2e cu)",
            e.hot[1], e.hot[2],
            e.cold[1], e.cold[2],
            e.Q, 
            e.ΔT, e.U, e.A, e.cost)
end
# Heat exchanger <<henheatexchanger>>:3 ends here

# [[file:../book.org::henprocesssteps][henprocesssteps]]
abstract type Step end

# the start of a hot/cold stream
struct Inlet <: Step
    inlet :: Stream
end

# the key step: exchange heat
# between hot and col streams
mutable struct Exchange <: Step
    hot :: Int
    cold :: Int
    match       # should be a node
    Q :: Float64
    hex :: Union{Exchanger, Nothing}
end

# no-operation: an empty step
struct Noop <: Step
end

# a set of steps, known as a segment,
# is a step itself and will be used
# for split & mix branches
mutable struct Segment <: Step
    steps :: Vector{Step}
end
# some segments are a single step
Segment(s :: Step) = Segment([s])

# split the stream
mutable struct Split <: Step
    η :: Union{Float64, Int}
end

# rejoin the split streams
struct Mix <: Step
end

# meet demand with utility
mutable struct Utility <: Step
    utility :: Union{ExternalUtility, Nothing}
    hex :: Union{Exchanger, Nothing}
end
# henprocesssteps ends here

# [[file:../book.org::hensegmentutilities][hensegmentutilities]]
function append!(segment :: Segment, step :: Step)
    push!(segment.steps, step)
end
function append!(segment :: Segment, steps :: Vector{Step})
    for s ∈ steps
        if ! (s isa Noop)
            append!(segment, s)
        end
    end
end
# hensegmentutilities ends here

# [[file:../book.org::henshowsteps][henshowsteps]]
function Base.show(io :: IO, segment :: Segment)
    print(io, " [")
    for s ∈ segment.steps
        print(io, s)
    end
    print(io, " ]")
end
Base.show(io :: IO, i :: Inlet) = print(io,"$(i.inlet)")
function Base.show(io :: IO, x :: Exchange)
    print(io, " x:$(x.hot)→$(x.cold)")
    if nothing != x.hex
        print(io, " $(x.hex)")
    else
        print(io, " Q=$(x.Q)")
    end
end
Base.show(io :: IO, s :: Split) = print(io, " ⊣ η=$(s.η)")
Base.show(io :: IO, m :: Mix) = print(io, " ⊢")
function Base.show(io :: IO, u :: Utility)
    if nothing != u.hex
        print(io, " u:$(u.hex)")
    else
        print(io, " u")
    end
end
Base.show(io :: IO, n :: Noop) = print(io, "")
# henshowsteps ends here

# [[file:../book.org::henparse][henparse]]
function parse(expression :: String, id :: Int)
    if nothing != (m = match(r"^ *(<|[0-9]+|\||>)",
                             expression))
        s = m.captures[1]
        expression = replace(expression,
                             r"^"*m.match => "")
        if s == "<"
            b1, expression = parse(expression, id)
            b2, expression = parse(expression, id)
            rest, expression = parse(expression, id)
            ( [ Split(0.0)   # split
                Segment(b1)  # branch 1
                Segment(b2)  # branch 2
                Mix()        # come back together
                rest ],      # remaining path
              expression )     
        elseif s == "|"
            # end of first branch
            ( Noop(), expression )
        elseif s == ">"
            # end of second branch
            ( Noop(), expression )
        else
            # must be a stream id
            n = Meta.parse(s)
            # assume hot streams have lower index
            # than cold streams
            from, to = id < n ? (id, n) : (n, id)
            rest, expression = parse(expression, id)
            ( vcat([Exchange(from, to, nothing,
                             0.0, nothing)],
                   rest),
              expression )
        end
    else
        ( Utility(nothing, nothing), expression )
    end
end
# henparse ends here

# [[file:../book.org::*Implementation of parser][Implementation of parser:2]]
function parse(expressions :: Vector{String})
    network = Segment[]
    for i ∈ 1:length(expressions)
        push!(network,
              Segment(parse(expressions[i], i)[1]))
    end
    network
end
# Implementation of parser:2 ends here

# [[file:../book.org::hennetworkgraphedge][hennetworkgraphedge]]
mutable struct Edge
    from       # node
    to         # node
    stream     # state
    id         # unique id
    traversed  # for evaluation
end
# hennetworkgraphedge ends here

# [[file:../book.org::hennetworkgraphnode][hennetworkgraphnode]]
mutable struct Node
    in         # edges
    out        # edges
    step       # task
    type       # :hot or :cold
    id         # unique id
    evaluated  # for costing
    z          # objective
    g          # infeasibility
end
# hennetworkgraphnode ends here

# [[file:../book.org::hennetworkgraph][hennetworkgraph]]
mutable struct Network
    inlets :: Vector{Node}  # root nodes
    outlets :: Vector{Node} # leaf nodes
    edges :: Vector{Edge}   # all edges in graph
    nodes :: Vector{Node}   # all nodes in graph
    ne :: Int               # next edge ID
    nn :: Int               # next node ID
    # empty network constructor
    Network() = new(Node[], Node[],
                    Edge[], Node[], 1, 1)
end
# hennetworkgraph ends here

# [[file:../book.org::hennetworkgraphconstructors][hennetworkgraphconstructors]]
function Edge(n :: Network)
    e = Edge(nothing, nothing, nothing, n.ne, false);
    add(n, e);
    e
end
function Edge(n :: Network, s :: Stream)
    e = Edge(nothing, nothing, s, n.ne, false);
    add(n, e);
    e
end
function Node(n :: Network, type)
    node = Node(nothing, nothing, nothing,
                type, n.nn,
                false,
                zeros(3),
                0.0);
    add(n, node);
    node
end
function Node(n :: Network, s :: Step, type)
    node = Node(nothing, nothing, s,
                type, n.nn,
                false,
                zeros(3),
                0.0);
    add(n, node);
    node
end
iscold(n :: Node) = n.type == :cold
ishot(n :: Node) = n.type == :hot
function add(network :: Network, edge :: Edge) 
    push!(network.edges, edge)
    network.ne += 1
end
function add(network :: Network, node :: Node) 
    push!(network.nodes, node)
    network.nn += 1
end
# hennetworkgraphconstructors ends here

# [[file:../book.org::segment2graph][segment2graph]]
function segment2graph(network :: Network,
                       edge :: Edge,
                       segment :: Segment,
                       type)
    s = 1                   # count through steps
    ns = length(segment.steps)
    first = true
    let root = nothing, leaf = nothing
        while s ≤ ns
            step = segment.steps[s]
            s += 1
            # ignore Noop steps
            if ! (step isa Noop)
                node = Node(network, step, type)
                edge.to = node
                if first
                    first = false
                    root = node
                end
                node.in = [edge]
                if step isa Exchange
                    edge = Edge(network)
                    node.out = [edge]
                    edge.from = node
                elseif step isa Split
                    # first branch
                    e1 = Edge(network)
                    e1.from = node
                    # next step should be a segment
                    step = segment.steps[s]
                    s += 1
                    # create a sub-graph from segment
                    (root, leaf1) = segment2graph(network, e1, step, type)
                    e1.to = root
                    # second branch
                    e2 = Edge(network)
                    e2.from = node
                    # connect branches to split node
                    node.out = [e1, e2]
                    step = segment.steps[s]
                    s += 1
                    # create a sub-graph from segment
                    (root, leaf2) = segment2graph(network, e2, step, type)
                    e2.to = root
                    # now join branches up: Mix
                    step = segment.steps[s]
                    s += 1
                    @assert step isa Mix "Expected a Mix after segments but got $step" 
                    node = Node(network, step, type)
                    # out edges of each branch leaf
                    # node are in-edges of the Mix
                    # node
                    node.in = vcat(leaf1.out,
                                   leaf2.out)
                    leaf1.out[1].to = node
                    leaf2.out[1].to = node
                    # prepare to connect to what
                    # follows
                    oldedge = edge
                    edge = Edge(network)
                    node.out = [edge]
                    edge.from = node
                elseif step isa Utility 
                    edge = Edge(network)
                    node.out = [edge]
                    edge.from = node
                else
                    error("We should not get here")
                end
                # the last node will be a leaf of
                # the directed graph
                leaf = node
            else
            end
        end
        # return root and leaf of segment
        return (root, leaf)
    end
end
# segment2graph ends here

# [[file:../book.org::hennetwork][hennetwork]]
function createnetwork(streams, segments)
    @assert length(streams) == length(segments) "Need a segment for each stream and vice versa."
    network = Network()
    for i ∈ 1:length(streams)
        # create root (inlet) node for stream
        node = Node(network,
                    Inlet(streams[i]),
                    streams[i].type)
        push!(network.inlets, node)
        # first edge from root node is the stream
        edge = Edge(network, streams[i])
        edge.from = node
        node.out = [edge]
        # and this edge is connected to the graph
        # corresponding to the segment for this
        # stream
        (root, leaf) = segment2graph(network, edge, segments[i], streams[i].type)
        # the leaf is also remembered to allow for
        # backward traversal
        push!(network.outlets, leaf)
    end
    # link up all process stream exchanges
    matchexchanges(network)
end
# hennetwork ends here

# [[file:../book.org::henmatch][henmatch]]
function matchexchanges(network)
    for node ∈ network.nodes
        if ishot(node) && node.step isa Exchange
            # for each hot exchange node, look for
            # corresponding cold node, in reverse
            # traversal of cold stream
            for i ∈ length(network.nodes):-1:1
                n = network.nodes[i]
                if iscold(n) && n.step isa Exchange
                    # for cold side, link up the two
                    # nodes if this cold node is not
                    # already linked
                    if n.step.hot == node.step.hot && n.step.cold == node.step.cold && nothing == n.step.match
                        # link up
                        n.step.match = node
                        node.step.match = n
                        break
                    end
                end
            end
            @assert node.step.match != nothing "Did not find a match"
        end
    end
    network
end
# henmatch ends here

# [[file:../book.org::hensproblem][hensproblem]]
mutable struct HENS
    streams :: Vector{Stream}
    representation :: Vector{String}
    network :: Network
    hexmodel # f(A)
    tolerance :: Real
    utilities :: Vector{ExternalUtility}
    function HENS(s, r, h, t, u)
        network = createnetwork(s,parse(r))
        new(s, r, network, h, t, u)
    end
end
# hensproblem ends here

# [[file:../book.org::hendesignvariables][hendesignvariables]]
getVariable(step :: Exchange) = step.Q
getVariable(step :: Split) = step.η
getVariables(problem :: HENS) = getVariables(problem.network)
function getVariables(network :: Network)
    [getVariable(n.step)
     for n ∈ network.nodes
         if ((ishot(n) && n.step isa Exchange)
             || n.step isa Split)]
end

setVariable(step :: Exchange, Q) = step.Q = Q
setVariable(step :: Split, η) = step.η = η
setVariables(problem :: HENS, vars) = setVariables(problem.network, vars)
function setVariables(network :: Network, vars)
    map((n,v) -> setVariable(n.step,v),
        [n for n ∈ network.nodes
             if ((ishot(n) && n.step isa Exchange)
                 || n.step isa Split)],
        vars)
end
# hendesignvariables ends here

# [[file:../book.org::henevaluatenetwork][henevaluatenetwork]]
function evaluate(vars, problem)
    network = problem.network
    setVariables(network, vars)
    # reset network: flags, costs, designs
    map(e -> (e.traversed = false), network.edges)
    for n in network.nodes
        n.evaluated = false
        n.z = [0.0, 0.0, 0.0]
        n.g = 0.0
        if n.step isa Exchange || n.step isa Utility
            n.step.hex = nothing
        end
    end
    # evaluate whole network with breadth first
    # traversal: hot streams first
    while reduce(|,
                 [(ishot(n) && !n.evaluated)
                  for n ∈ network.nodes])
        for n ∈ network.inlets
            if ishot(n) && !n.evaluated
                evaluate(n, problem)     # recursive
            end
        end
    end
    # now cold streams
    while reduce(|,
                 [(!ishot(n) && !n.evaluated)
                  for n ∈ network.nodes])
        for n ∈ network.inlets
            if !ishot(n) && !n.evaluated
                evaluate(n, problem)     # recursive
            end
        end
    end
    # collect and add up objective function and
    # feasibility values
    (reduce(+, [n.z for n ∈ network.nodes]),
     reduce(+, [n.g for n ∈ network.nodes]))
end
# henevaluatenetwork ends here

# [[file:../book.org::henevaluatenode][henevaluatenode]]
function evaluate(node :: Node, problem :: HENS)
    # check if node is ready for evaluation
    if node.in == nothing ||
        reduce(&, [edge.traversed
                   for edge ∈ node.in])
        # node can be evaluated which means designing
        # the processing step associated with it.
        design(node.step, node, problem)
        for i ∈ 1:length(node.out)
            node.out[i].traversed = true
            # recursively traverse graph
            if nothing != node.out[i].to
                evaluate(node.out[i].to, problem)
            end
        end
        node.evaluated = true
    end
end
# henevaluatenode ends here

# [[file:../book.org::hendesigninlet][hendesigninlet]]
function design(step :: Inlet, node, problem)
    # copy stream to the first (& only) edge
    node.out[1].stream = Stream(step.inlet)
end
# hendesigninlet ends here

# [[file:../book.org::hendesignexchange][hendesignexchange]]
function design(step :: Exchange, node, problem)
    @assert length(node.in) == 1 "Exchange $type node expects one inlet"
    in = node.in[1].stream      # inlet stream

    # actual exchange duty: if hot, it is defined by
    # the design variable (step.Q) as a fraction of
    # the available duty (in.Q); if a cold stream,
    # the value should have been set already when
    # the hot side was processed.
    if ishot(node)
        Q = step.Q * in.Q
        # tell cold side the amount to (potentially)
        # exchange
        step.match.step.Q = Q
    else
        Q = step.Q
    end
    # define the stream for the edge leaving this
    # node as being the same as the stream coming
    # in.  This stream will be modified by the
    # exchange if an exchange actually takes place.
    # Very small values of a desired exchange are
    # ignored.
    node.out[1].stream = Stream(in)

    # only do any work if we have an actual amount
    # to exchange; note that a negative value
    # indicates that the hot stream wished to give
    # away more heat than the cold stream can
    # accept.  This will be caught in the design of
    # the exchanger.
    if abs(Q) > problem.tolerance
        if Q ≤ in.Q
            T = in.Tin + (Q/in.Q)*(in.Tout-in.Tin)
            # duty for stream is adjusted
            node.out[1].stream.Q = in.Q-Q
            # as is the inlet temperature
            node.out[1].stream.Tin = T
            # if cold stream, design actual exchanger
            if ! ishot(node)
                try
                    node.step.hex = Exchanger(
                        step.match.in[1].stream,
                        in,
                        Q,
                        problem.hexmodel)
                    # three criteria
                    node.z = [node.step.hex.cost,
                              node.step.hex.A,
                              0.0]
                catch e
                    # the design may fail (Q<0, ΔT<0)
                    if e isa Infeasibility
                        node.g = e.g
                    else
                        # something else went wrong
                        # so propagate error upwards
                        throw(e)
                    end
                end
            end
        else
            # infeasible exchange: violation is
            # difference between heat required and
            # heat available
            node.g = abs(Q-in.Q)
        end
    end
end
# hendesignexchange ends here

# [[file:../book.org::hendesignutility][hendesignutility]]
function design(step :: Utility, node, problem)
    @assert length(node.in) == 1 "Exchange $type node expects one inlet"
    # inlet stream only
    in = node.in[1].stream
    # design an actual exchange if the amount to
    # heat/cool is significant.  This caters for
    # rounding errors in the processing of the heat
    # duties in getting to this node.
    if in.Q > problem.tolerance
        try
            # find appropriate utility; assume
            # utilities specified in order of
            # increasing cost so first suitable one
            # will do
            found = false
            for u ∈ problem.utilities
                if ishot(node) && !ishot(u) &&
                    in.Tout > u.Tout
                    found = true
                    node.step.hex = Exchanger(
                        in, u, in.Q,
                        problem.hexmodel)
                    utilitycost = u.model(in.Q)
                    # three criteria
                    node.z = [
                        node.step.hex.cost + utilitycost,
                        node.step.hex.A,
                        in.Q ]
                    # leave loop as utility found
                    break
                elseif !ishot(node) && ishot(u) &&
                    in.Tout < u.Tout
                    found = true
                    node.step.hex = Exchanger(
                        u, in, in.Q,
                        problem.hexmodel)
                    utilitycost = u.model(in.Q)
                    # three criteria
                    node.z = [
                        node.step.hex.cost + utilitycost,
                        node.step.hex.A,
                        in.Q ]
                    # leave loop as utility found
                    break
                end
            end
            if !found
                # no suitable external utility has
                # been found so we throw an
                # exception
                throw(Infeasibility(in.Q))
            end
        catch e
            if e isa Infeasibility
                node.g = e.g
            else
                # something else went wrong
                # so propagate error upwards
                throw(e)
            end
        end
    elseif in.Q < -problem.tolerance
        node.g = abs(in.Q)
    end
end
# hendesignutility ends here

# [[file:../book.org::hendesignsplit][hendesignsplit]]
function design(step :: Split, node, problem)
    in = node.in[1].stream
    mcp = [step.η, (1-step.η)] * in.mcp
    for i ∈ 1:2
        s = Stream(in)
        s.mcp = mcp[i]
        s.Q = abs(s.mcp*(s.Tout-s.Tin))
        node.out[i].stream = s
    end
end
# hendesignsplit ends here

# [[file:../book.org::hendesignmix][hendesignmix]]
function design(step :: Mix, node, problem)
    in = [node.in[1].stream, node.in[2].stream]
    node.out[1].stream = Stream(in[1])
    mcp = in[1].mcp + in[2].mcp
    Q = in[1].Q + in[2].Q
    if Q > eps()
        if ishot(node)
            Tin = in[1].Tout + Q/mcp
        else
            Tin = in[1].Tout - Q/mcp
        end
        node.out[1].stream.mcp = mcp
        node.out[1].stream.Q = Q
        node.out[1].stream.Tin = Tin
    else
        node.out[1].stream.Tin = in[1].Tout
    end
end
# hendesignmix ends here

# [[file:../book.org::hendesigninlet][hendesigninlet]]
function design(s :: Nothing, node, problem)
end
# hendesigninlet ends here

# [[file:../book.org::henevaluatenoop][henevaluatenoop]]
function design(step :: Noop, node, problem)
    node.out[1].stream = Stream(node.in[1].stream)
end
# henevaluatenoop ends here

# [[file:../book.org::hencostobjective][hencostobjective]]
function cost(v, problem)
    z, g = evaluate(v, problem)
    (z[1], g)
end
# hencostobjective ends here

# [[file:../book.org::henareaobjective][henareaobjective]]
function area(v, problem)
    z, g = evaluate(v, problem)
    (z[2], g)
end
# henareaobjective ends here

# [[file:../book.org::henutilobjective][henutilobjective]]
function utility(v, problem)
    z, g = evaluate(v, problem)
    (z[3], g)
end
# henutilobjective ends here

# [[file:../book.org::hennetworkshow][hennetworkshow]]
Base.show(io :: IO, e :: Edge) = e.stream == nothing ? print(io, "e$(e.id)<ϕ>") : print(io, "e$(e.id)<$(e.stream)>")
Base.show(io :: IO, n :: Node) = n.step == nothing ?
    print(io, "n$(n.id)<ϕ>") :
    (@printf(io, "n%d[(%.2e,%.2e,%.2e),%.1f]",
             n.id,
             n.z[1], n.z[2], n.z[3],
             n.g),
     print(io, " <$(n.step)>"))
function dump(desc :: String, edge :: Edge)
    println("$edge: $desc")
    println(": from = $(edge.from)")
    println(": to = $(edge.to)")
end
function dump(desc :: String, node :: Node)
    println("$node: $desc")
    println(": in = $(node.in)")
    println(": out = $(node.out)")
    println(": step = $(node.step)")
end
function traverse(io :: IO, node :: Node)
    if node isa Node
        print(io, node)
        if node.out == nothing
            print(io, "□")
        else
            for edge ∈ node.out
                print(io, edge)
                if edge.to == nothing
                    print(io, "?")
                else
                    traverse(io, edge.to)
                end
            end
        end
    else
        print(io, "□")        
    end
end
function Base.show(io :: IO, network :: Network)
    for node ∈ network.inlets
        traverse(io, node)
        println(io)
    end
end
# hennetworkshow ends here

# [[file:../book.org::henprintnetwork][henprintnetwork]]
using Printf
function printNetwork(network :: Network)
    n = length(network.nodes)
    c = Matrix{String}(undef, n, n)
    for i ∈ 1:n                 # for each node
        for j ∈ 1:n
            c[i,j] = "   "
        end
        node = network.nodes[i]
        if nothing != node.in && length(node.in) > 0
            if nothing != node.out && length(node.out) > 0 && any(e.to != nothing for e in node.out)
                c[i,i] = " +-"
            else
                c[i,i] = " # "
            end
        else
            c[i,i] = " >-"
        end
        for e ∈ node.out
            if nothing != e.to
                out = e.to
                c[i,out.id] = @sprintf "%2d" e.id
                for j ∈ i+1:out.id-1
                    c[i,j] = "---"
                    c[j,out.id] = "|"
                end
            end
        end
    end
    # top line of matrix
    @printf "     +-"
    for i ∈ 1:n
        @printf "%3d" i
    end
    @printf "-+\n"
    for i ∈ 1:n
        @printf "%3d. | " i
        for j ∈ 1:n
            @printf "%-3s" c[i,j]
        end
        @printf " | %s\n" network.nodes[i]
    end
    # bottom line of matrix
    @printf "     +-"
    for i ∈ 1:n
        @printf "%3d" i
    end
    @printf "-+\n"
    # now detail on the nodes and edges
    println("Nodes:")
    for node ∈ network.nodes
        @printf "%3d. %s\n" node.id node
    end
    println("Edges:")
    for edge ∈ network.edges
        @printf "%3d. %s\n" edge.id edge
    end
end
# henprintnetwork ends here

# [[file:../book.org::hendraw][hendraw]]
function traverse(n :: Node)
    print("    n$(n.id) ")
    if n.step isa Exchange
        if ishot(n)
            @printf " [label=\"X: %.3f Q\" shape=box]\n" n.step.Q
        else
            if n.step.Q > eps()
                @printf " [label=\"X: Q=%.1f\" shape=box]\n" n.step.Q
            else
                @printf " [label=\"X: Q=%.1f\" shape=box fontcolor=grey]\n" n.step.Q
            end
        end
    elseif n.step isa Utility
        if n.step.hex != nothing
            @printf(" [label=\"Utility\\nQ=%.1f\\nT=%.1f\\nA=%.2f\" shape=box color=%s]\n",
                    n.in[1].stream.Q,
                    n.in[1].stream.Tout,
                    n.step.hex.A,
                    (ishot(n) ? "blue" : "red"))
        else
            @printf(" [label=\"Utility\\nQ=%.1f\\nT=%.1f\" shape=box color=%s fontcolor=grey]\n",
                    n.in[1].stream.Q,
                    n.in[1].stream.Tout,
                    (ishot(n) ? "blue" : "red"))
        end
    elseif n.step isa Inlet
        @printf " [label=\"Q=%.1f\\nT=%.1f\" shape=plaintext]\n" n.step.inlet.Q n.step.inlet.Tin
    elseif n.step isa Mix
        @printf " [label=\"mix\" shape=plaintext]\n"
    elseif n.step isa Split
        @printf(" [label=\"split: η=%.3f\\nmcp=%.3f/%.3f\" shape=plaintext]\n",
                n.step.η,
                n.step.η*n.in[1].stream.mcp,
                (1-n.step.η)*n.in[1].stream.mcp)
    else
        @printf " [label=\"\" shape=plaintext]\n"
    end
    # traverse this graph depth first
    if ishot(n)
        if n.out != nothing
            for e in n.out
                if e.to != nothing
                    traverse(e.to)
                end
            end
        end
    else
        if n.in != nothing
            for e in n.in
                if e.from != nothing
                    traverse(e.from)
                end
            end
        end
    end
end
function draw(problem :: HENS)
    network = problem.network
    # start of graphviz
    println("#+begin_src dot :file network.png :results file :eval no-export")
    println("graph network {
      overlap=scale
      nodesep=1
      ranksep=1 ")
    # reset edges: flags, costs, designs
    map(e -> (e.traversed = false), network.edges)
    s = 1                       # counter for streams
    # hot streams
    for n in network.inlets
        if ishot(n)
            println("  subgraph cluster$s { color=red; label=\"$(problem.streams[s])\"")
            s += 1
            traverse(n)
            println("  }")
        end
    end
    # cold streams now
    for n in network.outlets
        if ! ishot(n)
            println("  subgraph cluster$s { color=blue; label=\"$(problem.streams[s])\"")
            s += 1
            traverse(n)
            println("  }")
        end
    end
    # now connect the nodes, including exchanges
    for e in network.edges
        if e.from != nothing && e.to != nothing
            if ishot(e.from)
                print("  n$(e.from.id) -- n$(e.to.id)")
            else
                print("  n$(e.to.id) -- n$(e.from.id)")
            end
            if e.stream != nothing
                if e.stream.Q > eps()
                    @printf(" [label=\"Q=%.1f\\nT=%.1f\"]\n",
                            e.stream.Q, e.stream.Tin)
                else
                    @printf(" [label=\"Q=%.1f\\nT=%.1f\" fontcolor=grey]\n",
                            e.stream.Q, e.stream.Tin)
                end
            else
                println()
            end
        end
    end

    # finally, the exchanges between hot and cold
    for n in network.nodes
        if ishot(n)
            if n.step isa Exchange
                cold = n.step.match
                hex = cold.step.hex
                if hex != nothing
                    @printf("  n%d:e -- n%d:e
                            [label=\"A=%.2f\",
                            style=dashed,
                            constraint=false,
                            weight=0]\n",
                            n.id,
                            n.step.match.id,
                            hex.A)
                else
                    @printf("  n%d:e -- n%d:e 
                         [label=\"A=0\",
                         style=dashed,
                            constraint=false, 
                            weight=0,
                            color=grey,
                            fontcolor=grey]\n",
                            n.id,
                            n.step.match.id)
                end
            end
        end
    end
    println("}")
    println("#+end_src")
end
# hendraw ends here
