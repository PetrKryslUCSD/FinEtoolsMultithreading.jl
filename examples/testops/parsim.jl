println("Current folder: $(pwd())")

if length(ARGS) < 1
    error("I need at least one arguments: N (mesh subdivision)")
end

using Pkg

# Pkg.add("ThreadPinning")

Pkg.activate(".")
Pkg.instantiate()

using LinearAlgebra
LinearAlgebra.BLAS.set_num_threads(1)

# Turn off thread pinning because it seems to interfere with the graph coloring library.
# using ThreadPinning
# ThreadPinning.Prefs.set_os_warning(false)
# pinthreads(:cores)

N = parse(Int, ARGS[1])
ntasks = Threads.nthreads()
if length(ARGS) > 1
    ntasks = parse(Int, ARGS[2])
end
assembly_only = true
if length(ARGS) > 2
    assembly_only = parse(Bool, ARGS[3])
end

include(raw"ops_parallel.jl")
using .ops_parallel; 

using FinEtoolsMultithreading
Pkg.status("FinEtoolsMultithreading")

NTRIALS = 5
for trial in 1:NTRIALS
    @info "Trial $(trial) out of $(NTRIALS): nthreads=$(Threads.nthreads()), ntasks=$(ntasks), N=$(N)"
    ops_parallel.run(N, ntasks,  assembly_only)
    GC.gc(true)
end

