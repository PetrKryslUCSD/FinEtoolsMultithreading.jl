module stubby_corbel_serial
using FinEtools
using FinEtools.AlgoBaseModule: evalconvergencestudy, solve_blocked!
using FinEtoolsDeforLinear
using FinEtoolsDeforLinear.AlgoDeforLinearModule:
    linearstatics, exportstresselementwise, exportstress
using Statistics: mean
using LinearAlgebra
using SparseArrays
using SuiteSparse
using SymRCM
using Random        
using DataDrop

# Isotropic material
E = 1000.0
nu = 0.3 # Compressible material
W = 25.0
H = 50.0
L = 50.0
htol = minimum([L, H, W]) / 1000
uzex = -2.16e-01
magn = 0.2 * (-12.6) / 4
Force = magn * W * H * 2
CTE = 0.0

function getfrcL!(forceout, XYZ, tangents, feid, qpid)
    copyto!(forceout, [0.0; 0.0; magn])
end

function run(N = 10, assembly_only = false)
    times = Dict{String,Vector{Float64}}()
    
    t1 = time()
    fens, fes = H8block(W, L, H, N, 2 * N, 2 * N)
    times["MeshGeneration"] = [time() - t1]
    println("Number of elements: $(count(fes))")
    bfes = meshboundary(fes)
    # end cross-section surface  for the shear loading
    sectionL = selectelem(fens, bfes; facing = true, direction = [0.0 +1.0 0.0])
    # 0 cross-section surface  for the reactions
    section0 = selectelem(fens, bfes; facing = true, direction = [0.0 -1.0 0.0])
    # 0 cross-section surface  for the reactions
    sectionlateral = selectelem(fens, bfes; facing = true, direction = [1.0 0.0 0.0])

    MR = DeforModelRed3D
    material = MatDeforElastIso(MR, 0.0, E, nu, CTE)

    # Material orientation matrix
    csmat = [i == j ? one(Float64) : zero(Float64) for i = 1:3, j = 1:3]

    function updatecs!(csmatout, XYZ, tangents, feid, qpid)
        copyto!(csmatout, csmat)
    end

    geom = NodalField(fens.xyz)
    u = NodalField(zeros(size(fens.xyz, 1), 3)) # displacement field

    lx0 = connectednodes(subset(bfes, section0))
    setebc!(u, lx0, true, 1, 0.0)
    setebc!(u, lx0, true, 2, 0.0)
    setebc!(u, lx0, true, 3, 0.0)
    lx1 = connectednodes(subset(bfes, sectionlateral))
    setebc!(u, lx1, true, 1, 0.0)
    applyebc!(u)
    numberdofs!(u)
    # numberdofs!(u)
    println("nfreedofs(u) = $(nfreedofs(u))")

    fi = ForceIntensity(Float64, 3, getfrcL!)
    el2femm = FEMMBase(IntegDomain(subset(bfes, sectionL), GaussRule(2, 2)))
    F = distribloads(el2femm, geom, u, fi, 2)
    F_f = vector_blocked_f(F, nfreedofs(u))

    # femm = FEMMDeforLinear(MR, IntegDomain(fes, GaussRule(3, 2)), material)
    femm = FEMMDeforLinearMSH8(MR, IntegDomain(fes, GaussRule(3, 2)), material)

    t0 = time()
    t1 = time()
    associategeometry!(femm, geom)
    assembler = SysmatAssemblerSparse(1.0)  
    setnomatrixresult(assembler, true)
    stiffness(femm, assembler, geom, u)
    times["ComputeCOO"] = [time() - t1]
    println("Compute COO = $(times["ComputeCOO"]) [s]")
    t1 = time()
    setnomatrixresult(assembler, false)
    K = makematrix!(assembler)
    times["BuildCSC"] = [time() - t1]
    println("Build CSC = $(times["BuildCSC"]) [s]")

    times["TotalAssembly"] = [time() - t0]
    println("Assembly total = $(times["TotalAssembly"]) [s]")

    K_ff = matrix_blocked_ff(K, nfreedofs(u))
    F_f = vector_blocked_f(F, nfreedofs(u))
    println("Stiffness: number of non zeros = $(nnz(K_ff)) [ND]")
    println("Sparsity = $(nnz(K_ff)/size(K_ff, 1)/size(K_ff, 2)) [ND]")
    
    if assembly_only
        isdir("$(N)") || mkdir("$(N)")
        n = DataDrop.with_extension(joinpath("$(N)", "stubby_corbel_H8_ms-timing-serial"), "json")
        if isfile(n)
            storedtimes = DataDrop.retrieve_json(n)
            for k in keys(storedtimes)
                times[k] = cat(times[k], storedtimes[k], dims = 1)
            end
        end
        DataDrop.store_json(n, times)
        return
    end
    
    Tipl = selectnode(fens, box = [0 W L L 0 H], inflate = htol)
    @show norm(K_ff - K_ff') / norm(K_ff)
    @time K_ff_factors = SuiteSparse.CHOLMOD.cholesky(Symmetric(K_ff))
    @show nnz(K_ff_factors)
    @time U_f = K_ff_factors \ (F_f)

    scattersysvec!(u, U_f[:])

    utip = mean(u.values[Tipl, 3], dims = 1)
    println("Deflection: $(utip), compared to $(uzex)")

    File = "stubby_corbel_H8_ms_serial.vtk"
    vtkexportmesh(File, fens, fes; vectors = [("u", u.values)])
    
    true
end # run

end # module 
nothing
