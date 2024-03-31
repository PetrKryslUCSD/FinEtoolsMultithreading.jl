module sphere_modes_parallel
using FinEtools
using FinEtools.AlgoBaseModule: matrix_blocked
using FinEtoolsAcoustics
using FinEtoolsMultithreading
using FinEtoolsMultithreading.Exports
using FinEtoolsMultithreading: domain_decomposition, 
          parallel_matrix_assembly!, SysmatAssemblerSparsePatt, SysmatAssemblerSparsePattwLookup
using FinEtools.MeshExportModule
using LinearAlgebra
using Arpack: eigs
using DataDrop

# For the data 
# rho = 1.2*phun("kg/m^3");# mass density
# c  = 340.0*phun("m/s");# sound speed
# bulk =  c^2*rho;
# R = 1000.0*phun("mm");# radius of the piston

# the reference 

# @article{GAO2013914,
# title = {Eigenvalue analysis for acoustic problem in 3D by boundary element method with the block Sakurai–Sugiura method},
# journal = {Engineering Analysis with Boundary Elements},
# volume = {37},
# number = {6},
# pages = {914-923},
# year = {2013},
# issn = {0955-7997},
# doi = {https://doi.org/10.1016/j.enganabound.2013.03.015},
# url = {https://www.sciencedirect.com/science/article/pii/S0955799713000714},
# author = {Haifeng Gao and Toshiro Matsumoto and Toru Takahashi and Hiroshi Isakari},
# keywords = {Eigenvalues, Acoustic, The block SS method, Boundary element method, Burton–Miller's method},
# abstract = {This paper presents accurate numerical solutions for nonlinear eigenvalue analysis of three-dimensional acoustic cavities by boundary element method (BEM). To solve the nonlinear eigenvalue problem (NEP) formulated by BEM, we employ a contour integral method, called block Sakurai–Sugiura (SS) method, by which the NEP is converted to a standard linear eigenvalue problem and the dimension of eigenspace is reduced. The block version adopted in present work can also extract eigenvalues whose multiplicity is larger than one, but for the complex connected region which includes a internal closed boundary, the methodology yields fictitious eigenvalues. The application of the technique is demonstrated through the eigenvalue calculation of sphere with unique homogenous boundary conditions, cube with mixed boundary conditions and a complex connected region formed by cubic boundary and spherical boundary, however, the fictitious eigenvalues can be identified by Burton–Miller's method. These numerical results are supported by appropriate convergence study and comparisons with close form.}
# }

# shows the wave numbers in Table 1.
#=
The multiplicity of the Dirichlet eigenvalues.
Wavenumber*R                    Multiplicity
3.14159, 6.28319, 9.42478       1
4.49340, 7.72525, 10.90412      3
5.76346, 9.09501, 12.32294      5
6.98793, 10.41711, 13.69802     7
8.18256, 11.70491, 15.03966     9
9.35581, 12.96653, 16.35471    11
=#

# Sphere with Dirichlet boundary conditions: model analysis.
# Sphere of radius $(R), in WATER.
# Tetrahedral T4 mesh.
# Exact fundamental frequency: $(c/2/R)
function run(N=2, ntasks=Threads.nthreads(), assembly_only=false)
    

    rho = 1000 * phun("kg/m^3")# mass density
    c = 1500.0 * phun("m/s")# sound speed
    bulk = c^2 * rho
    R = 500.0 * phun("mm")# radius of the sphere
    tolerance = R / 1e3

    neigvs = 7

    wn_table = [
        ([3.14159, 6.28319, 9.42478], 1),
        ([4.49340, 7.72525, 10.90412], 3),
        ([5.76346, 9.09501, 12.32294], 5),
        ([6.98793, 10.41711, 13.69802], 7),
        ([8.18256, 11.70491, 15.03966], 9),
        ([9.35581, 12.96653, 16.35471], 11),
    ]
    # @info "Reference frequencies"
    # for i in axes(wn_table, 1)
    #     fq = wn_table[i][1] ./ R .* c / (2 * pi)
    #     @info "$(fq), multiplicity $(wn_table[i][2])"
    # end

    fens, fes = H8sphere(R, N)
    renumb(c) = c[[1, 4, 3, 2, 5, 8, 7, 6]]
    fens1, fes1 = mirrormesh(
        fens,
        fes,
        [-1.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
        renumb=renumb,
    )
    fens, newfes1, fes2 = mergemeshes(fens1, fes1, fens, fes, tolerance)
    fes = cat(newfes1, fes2)
    fens1, fes1 = mirrormesh(
        fens,
        fes,
        [0.0, -1.0, 0.0],
        [0.0, 0.0, 0.0],
        renumb=renumb,
    )
    fens, newfes1, fes2 = mergemeshes(fens1, fes1, fens, fes, tolerance)
    fes = cat(newfes1, fes2)
    fens1, fes1 = mirrormesh(
        fens,
        fes,
        [0.0, 0.0, -1.0],
        [0.0, 0.0, 0.0],
        renumb=renumb,
    )
    fens, newfes1, fes2 = mergemeshes(fens1, fes1, fens, fes, tolerance)
    fes = cat(newfes1, fes2)

    @info "$(count(fens)) nodes"

    geom = NodalField(fens.xyz)
    P = NodalField(zeros(size(fens.xyz, 1), 1))
    bfes = meshboundary(fes)
    setebc!(P, connectednodes(bfes))
    numberdofs!(P)

    mass_times = Dict{String,Vector{Float64}}()

    t1 = time()
    n2e = FENodeToFEMap(fes.conn, nnodes(P))
    mass_times["FENodeToFEMap"] = [time() - t1]
    println("Make node to element map = $(mass_times["FENodeToFEMap"]) [s]")

    material = MatAcoustFluid(bulk, rho)

    GC.enable(false)

    t0 = time()

    t1 = time()
    e2e = FEElemToNeighborsMap(n2e, fes)
    mass_times["FEElemToNeighborsMap"] = [time() - t1]
    println("    Make element to neighbor map = $(mass_times["FEElemToNeighborsMap"]) [s]")

    t1 = time()
    coloring = FinEtoolsMultithreading.element_coloring(fes, e2e)
    mass_times["ElementColors"] = [time() - t1]
    println("    Compute element colors = $(mass_times["ElementColors"]) [s]")

    t1 = time()
    n2n = FENodeToNeighborsMap(n2e, fes)
    mass_times["FENodeToNeighborsMap"] = [time() - t1]
    println("    Make node to neighbor map = $(mass_times["FENodeToNeighborsMap"]) [s]")

    t1 = time()
    K_pattern = sparse_symmetric_csc_pattern(P.dofnums, nalldofs(P), n2n, zero(eltype(P.values)))
    mass_times["SparsityPattern"] = [time() - t1]
    println("    Sparsity pattern = $(mass_times["SparsityPattern"]) [s]")

    t1 = time()
    decomposition = domain_decomposition(fes, coloring,
        (fessubset) -> FEMMAcoust(IntegDomain(fessubset, GaussRule(3, 2)), material), ntasks)
    mass_times["DomainDecomposition"] = [time() - t1]
    println("    Domain decomposition = $(mass_times["DomainDecomposition"]) [s]")

    AT = SysmatAssemblerSparsePatt
    AT = SysmatAssemblerSparsePattwLookup

    @time AT(K_pattern)

    t1 = time()
    Ma = parallel_matrix_assembly!(
        AT(K_pattern),
        decomposition,
        (femm, assmblr) -> acousticmass(femm, assmblr, geom, P),
    )
    mass_times["AssemblyOfValues"] = [time() - t1]
    println("    Add to matrix = $(mass_times["AssemblyOfValues"]) [s]")

    mass_times["TotalAssembly"] = [time() - t0]
    println("Assembly MASS total = $(mass_times["TotalAssembly"]) [s]")

    GC.enable(true)

    GC.enable(false)

    stiffness_times = Dict{String,Vector{Float64}}()
    
    t0 = time()

    t1 = time()
    e2e = FEElemToNeighborsMap(n2e, fes)
    stiffness_times["FEElemToNeighborsMap"] = [time() - t1]
    println("    Make element to neighbor map = $(stiffness_times["FEElemToNeighborsMap"]) [s]")

    t1 = time()
    coloring = FinEtoolsMultithreading.element_coloring(fes, e2e)
    stiffness_times["ElementColors"] = [time() - t1]
    println("    Compute element colors = $(stiffness_times["ElementColors"]) [s]")

    t1 = time()
    n2n = FENodeToNeighborsMap(n2e, fes)
    stiffness_times["FENodeToNeighborsMap"] = [time() - t1]
    println("    Make node to neighbor map = $(stiffness_times["FENodeToNeighborsMap"]) [s]")

    t1 = time()
    K_pattern = sparse_symmetric_csc_pattern(P.dofnums, nalldofs(P), n2n, zero(eltype(P.values)))
    stiffness_times["SparsityPattern"] = [time() - t1]
    println("    Sparsity pattern = $(stiffness_times["SparsityPattern"]) [s]")

    t1 = time()
    decomposition = domain_decomposition(fes, coloring,
        (fessubset) -> FEMMAcoust(IntegDomain(fessubset, GaussRule(3, 2)), material), ntasks)
    stiffness_times["DomainDecomposition"] = [time() - t1]
    println("    Domain decomposition = $(stiffness_times["DomainDecomposition"]) [s]")

    t1 = time()
    Ka = parallel_matrix_assembly!(
        AT(K_pattern),
        decomposition,
        (femm, assmblr) -> acousticstiffness(femm, assmblr, geom, P),
    )
    stiffness_times["AssemblyOfValues"] = [time() - t1]
    println("    Add to matrix = $(stiffness_times["AssemblyOfValues"]) [s]")

    stiffness_times["TotalAssembly"] = [time() - t0]
    println("Assembly STIFFNESS total = $(stiffness_times["TotalAssembly"]) [s]")

    GC.enable(true)
 
    if assembly_only
        isdir("$(N)") || mkdir("$(N)")
        n = DataDrop.with_extension(joinpath("$(N)", "sphere_modes_parallel-timing-parallel-stiffness-nth=$(ntasks)"), "json")
        if isfile(n)
            storedtimes = DataDrop.retrieve_json(n)
            for k in keys(storedtimes)
                stiffness_times[k] = cat(stiffness_times[k], storedtimes[k], dims=1)
            end
        end
        DataDrop.store_json(n, stiffness_times)
        n = DataDrop.with_extension(joinpath("$(N)", "sphere_modes_parallel-timing-parallel-mass-nth=$(ntasks)"), "json")
        if isfile(n)
            storedtimes = DataDrop.retrieve_json(n)
            for k in keys(storedtimes)
                mass_times[k] = cat(mass_times[k], storedtimes[k], dims=1)
            end
        end
        DataDrop.store_json(n, mass_times)
        return
    end

    Ma_ff = matrix_blocked(Ma, nfreedofs(P), nfreedofs(P))[:ff]
    Ka_ff = matrix_blocked(Ka, nfreedofs(P), nfreedofs(P))[:ff]

    d, v, nconv = eigs(Ka_ff, Ma_ff; nev=neigvs, which=:SM, explicittransform=:none)
    v = real.(v)
    fs = real(sqrt.(complex(d))) ./ (2 * pi)
    @info("Frequencies (1:5): $(fs[1:5]) [Hz]")
    @info "Reference frequencies"
    for i in axes(wn_table, 1)
        fq = wn_table[i][1] ./ R .* c / (2 * pi)
        @info "$(fq), multiplicity $(wn_table[i][2])"
    end
    ks = (2 * pi) .* fs ./ c ./ phun("m")
    # @info("Wavenumbers: $(ks) [m]")

    File = "sphere_modes_parallel.vtk"
    scalarllist = Any[]
    for n = [2, 5, 7]
        scattersysvec!(P, v[:, n])
        push!(scalarllist, ("Pressure_mode_$n", deepcopy(P.values)))
    end
    vtkexportmesh(
        File,
        connasarray(fes),
        geom.values,
        FinEtools.MeshExportModule.VTK.H8;
        scalars=scalarllist,
    )
    # @async run(`"paraview.exe" $File`)

    true

end # sphere_h8_in_air

end # module sphere_mode_examples
nothing
