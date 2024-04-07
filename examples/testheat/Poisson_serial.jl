module Poisson_serial
using FinEtools
using FinEtools.AlgoBaseModule: solve_blocked!, matrix_blocked, vector_blocked
using FinEtools.AssemblyModule
using FinEtools.MeshExportModule
using FinEtoolsHeatDiff
using ChunkSplitters
using LinearAlgebra
using DataDrop
using SparseArrays
using SymRCM

const A = 1.0 # dimension of the domain (length of the side of the square)

function make_model_H20(N = 25)
    fens, fes = H20block(A, A, A, N, N, N)
    ir = GaussRule(3, 3)
    return fens, fes, ir
end

function _run(make_model, N, assembly_only)
    times = Dict{String,Vector{Float64}}()
    thermal_conductivity = [i == j ? one(Float64) : zero(Float64) for i = 1:3, j = 1:3] # conductivity matrix
    Q = -6.0 # internal heat generation rate
    function getsource!(forceout, XYZ, tangents, feid, qpid)
        forceout[1] = Q #heat source
    end
    tempf(x) = (1.0 .+ x[:, 1] .^ 2 + 2.0 .* x[:, 2] .^ 2)#the exact distribution of temperature1
    fens, fes, ir = make_model(N)
    geom = NodalField(fens.xyz)
    Temp = NodalField(zeros(size(fens.xyz, 1), 1))
    Tolerance = 1.0 / count(fes) / 100
    l1 = selectnode(fens; box = [0.0 0.0 0.0 A 0.0 A], inflate = Tolerance)
    l2 = selectnode(fens; box = [A A 0.0 A 0.0 A], inflate = Tolerance)
    l3 = selectnode(fens; box = [0.0 A 0.0 0.0 0.0 A], inflate = Tolerance)
    l4 = selectnode(fens; box = [0.0 A A A 0.0 A], inflate = Tolerance)
    l5 = selectnode(fens; box = [0.0 A 0.0 A 0.0 0.0], inflate = Tolerance)
    l6 = selectnode(fens; box = [0.0 A 0.0 A A A], inflate = Tolerance)
    List = vcat(l1, l2, l3, l4, l5, l6)
    setebc!(Temp, List, true, 1, tempf(geom.values[List, :])[:])
    numberdofs!(Temp)
    println("Number of free degrees of freedom: $(nfreedofs(Temp))")
    
    material = MatHeatDiff(thermal_conductivity)
    femm = FEMMHeatDiff(IntegDomain(fes, ir), material)
    println("Conductivity")
    t0 = time()
    t1 = time()
    K = conductivity(femm, geom, Temp)
    times["TotalAssemblyStiffness"] = [time() - t0]
    println("Assembly total = $(times["TotalAssemblyStiffness"]) [s]")

    if assembly_only
        isdir("$(N)") || mkdir("$(N)")
        n = DataDrop.with_extension(joinpath("$(N)", "Poisson_serial-timing-stiffness"), "json")
        if isfile(n)
            storedtimes = DataDrop.retrieve_json(n)
            for k in keys(storedtimes)
                times[k] = cat(times[k], storedtimes[k], dims=1)
            end
        end
        return
    end

    println("Internal heat generation")
    fi = ForceIntensity(Float64[Q])
    F1 = distribloads(femm, geom, Temp, fi, 3)
    println("Solution of the system")
    t1 = time()
    solve_blocked!(Temp, K, F1)
    times["Solution"] = [time() - t0]
    println("Solution = $(times["Solution"]) [s]")

    Error = 0.0
    for k in axes(fens.xyz, 1)
        Error = Error + abs.(Temp.values[k, 1] - tempf(reshape(fens.xyz[k, :], (1, 3)))[1])
    end
    println("Error =$Error")
    true
end # Poisson_serial

function run(N = 25, assembly_only = false)
    return _run(make_model_H20, N, assembly_only)
end

end # module Poisson_serial
nothing
