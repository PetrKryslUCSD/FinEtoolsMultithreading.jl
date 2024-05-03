module sphere_modes_parallel
using FinEtools
using FinEtools.AlgoBaseModule: matrix_blocked
using FinEtoolsAcoustics
using FinEtoolsMultithreading.Exports
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
function run(N = 2, ntasks = Threads.nthreads(), assembly_only = false)
    times = Dict{String, Vector{Float64}}()
    
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
        renumb = renumb,
    )
    fens, newfes1, fes2 = mergemeshes(fens1, fes1, fens, fes, tolerance)
    fes = cat(newfes1, fes2)
    fens1, fes1 = mirrormesh(
        fens,
        fes,
        [0.0, -1.0, 0.0],
        [0.0, 0.0, 0.0],
        renumb = renumb,
    )
    fens, newfes1, fes2 = mergemeshes(fens1, fes1, fens, fes, tolerance)
    fes = cat(newfes1, fes2)
    fens1, fes1 = mirrormesh(
        fens,
        fes,
        [0.0, 0.0, -1.0],
        [0.0, 0.0, 0.0],
        renumb = renumb,
    )
    fens, newfes1, fes2 = mergemeshes(fens1, fes1, fens, fes, tolerance)
    fes = cat(newfes1, fes2)

    @info "$(count(fens)) nodes"

    geom = NodalField(fens.xyz)
    P = NodalField(zeros(size(fens.xyz, 1), 1))
    bfes = meshboundary(fes)
    setebc!(P, connectednodes(bfes))
    numberdofs!(P)

    t1 = time()
    n2e = FENodeToFEMapThr(fes, nnodes(P))
    times["FENodeToFEMap"] = [time() - t1]
    println("Make node to element map = $(times["FENodeToFEMap"]) [s]")

   material =  MatAcoustFluid(bulk, rho)

   t1 = time()
    Ma = parallel_make_matrix(
        fes,
        P.dofnums,
        nalldofs(P),
        eltype(P.values),
        n2e,
        (fessubset) -> FEMMAcoust(IntegDomain(fessubset, GaussRule(3, 2)), material),
        (femm, assmblr) -> acousticmass(femm, assmblr, geom, P),
        ntasks,
        :CSC
    )
    # Ma = acousticmass(femm, geom, P)
    times["AssembleMass"] = [time() - t1]
    println("Assemble mass = $(times["AssembleMass"]) [s]")

    t1 = time()
    Ka = parallel_make_matrix(
        fes,
        P.dofnums,
        nalldofs(P),
        eltype(P.values),
        n2e,
        (fessubset) -> FEMMAcoust(IntegDomain(fessubset, GaussRule(3, 2)), material),
        (femm, assmblr) -> acousticstiffness(femm, assmblr, geom, P),
        ntasks,
        :CSC
    )
    # Ka = acousticstiffness(femm, geom, P)
    times["AssembleStiffness"] = [time() - t1]
    println("Assemble stiffness = $(times["AssembleStiffness"]) [s]")

    if assembly_only
        isdir("$(N)") || mkdir("$(N)")
        n = DataDrop.with_extension(joinpath("$(N)", "sphere_modes_parallel-timing-parallel-nth=$(ntasks)"), "json")
        if isfile(n)
            storedtimes = DataDrop.retrieve_json(n)
            for k in keys(storedtimes)
                times[k] = cat(times[k], storedtimes[k], dims = 1)
            end
        end
        DataDrop.store_json(n, times)
        return
    end

    Ma_ff = matrix_blocked(Ma, nfreedofs(P), nfreedofs(P))[:ff]
    Ka_ff = matrix_blocked(Ka, nfreedofs(P), nfreedofs(P))[:ff]

    d, v, nconv = eigs(Ka_ff, Ma_ff; nev = neigvs, which = :SM, explicittransform = :none)
    v = real.(v)
    fs = real(sqrt.(complex(d))) ./ (2 * pi)
    # @info("Frequencies (1:5): $(fs[1:5]) [Hz]")
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
        scalars = scalarllist,
    )
    # @async run(`"paraview.exe" $File`)

    true

end # sphere_h8_in_air

end # module sphere_mode_examples
nothing
