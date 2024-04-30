println("Current folder: $(pwd())")

if length(ARGS) < 1
    error("I need at least one arguments: N (mesh subdivision)")
end

using Pkg

Pkg.add("ThreadPinning")

Pkg.activate(".")
Pkg.instantiate()

using LinearAlgebra
LinearAlgebra.BLAS.set_num_threads(1)


using ThreadPinning
ThreadPinning.Prefs.set_os_warning(false)
pinthreads(:cores)

N = parse(Int, ARGS[1])
ntasks = Threads.nthreads()
if length(ARGS) > 1
    ntasks = parse(Int, ARGS[2])
end
assembly_only = true
if length(ARGS) > 2
    assembly_only = parse(Bool, ARGS[3])
end

include(raw"sphere_modes_parallel.jl")
using .sphere_modes_parallel; 

NTRIALS = 5
for trial in 1:NTRIALS
    @info "Trial $(trial) out of $(NTRIALS): nthreads=$(Threads.nthreads()), ntasks=$(ntasks), N=$(N)"
    sphere_modes_parallel.run(N, ntasks,  assembly_only)
    GC.gc(true)
end

