println("Current folder: $(pwd())")

if length(ARGS) < 1
    error("I need one argument: N (mesh subdivision)")
end

using Pkg

Pkg.activate(".")
Pkg.instantiate()

N = parse(Int, ARGS[1])
assembly_only = true
if length(ARGS) > 1
    assembly_only = parse(Bool, ARGS[2])
end

include(raw"sphere_modes_serial.jl")
using .sphere_modes_serial; 

NTRIALS = 5
for trial in 1:NTRIALS
    @info "Trial $(trial) out of $(NTRIALS): N=$(N)"
    sphere_modes_serial.run(N, assembly_only)
    GC.gc(true)
end

