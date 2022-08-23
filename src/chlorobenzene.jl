# [[file:../book.org::jacarandamodule][jacarandamodule]]
using Jacaranda
# jacarandamodule ends here

# [[file:../book.org::chlorobenzeneprocess][chlorobenzeneprocess]]
function process()
    Jacaranda.init("clben-process-flowsheet.in", "clben")
end
# chlorobenzeneprocess ends here

# [[file:../book.org::chlorobenzeneevaluation][chlorobenzeneevaluation]]
function evaluate(v, flowsheet)
    f = Jacaranda.evaluate(flowsheet, v)
    # Jacaranda returns, for this process, two
    # objective function values and three
    # constraint violations so put these in the
    # form that an optimization system can use
    (f[1:2], maximum(f[3:5]))
end
# chlorobenzeneevaluation ends here
