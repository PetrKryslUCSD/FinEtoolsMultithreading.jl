module mparallelassembly_assembler_1
using FinEtools
using FinEtoolsMultithreading.Exports
using FinEtoolsMultithreading: SysmatAssemblerSparsePatt, startassembly!, assemble!, makematrix!, associate_pattern
using LinearAlgebra
using Test

function test()
    W = 1.1
    L = 12.0
    t = 0.32
    nl, nt, nw = 12, 33, 24
    nl, nt, nw = 2, 3, 2
    ntasks = 2

    fens, fes = H8block(L, W, t, nl, nw, nt)
    geom = NodalField(fens.xyz)
    psi = NodalField(fill(1.0, count(fens), 1))
    nl = collect(1:3)
    setebc!(psi, nl, true, ones(Int, length(nl)), 0.0)
    numberdofs!(psi)
    v_f = gathersysvec(psi)

    femm = FEMMBase(IntegDomain(fes, GaussRule(3, 2)))
    n2e = FENodeToFEMap(fes.conn, count(fens))
    n2n = FENodeToNeighborsMap(n2e, fes.conn)
    K = sparse_symmetric_csc_pattern(psi.dofnums, nalldofs(psi), n2n, zero(eltype(psi.values)))
    assembler = SysmatAssemblerSparsePatt(0.0)
    associate_pattern(assembler, K)
    K = bilform_diffusion(femm, assembler, geom, psi, DataCache(Matrix(1.0 * LinearAlgebra.I(3))))
    
    K_ff = matrix_blocked_ff(K, nfreedofs(psi))
    result = abs(v_f' * K_ff * v_f)

    ass = SysmatAssemblerFFBlock(nfreedofs(psi))
    K_ff2 = bilform_diffusion(FEMMBase(IntegDomain(fes, GaussRule(3, 2))), ass, geom, psi, DataCache(Matrix(1.0 * LinearAlgebra.I(3))))
    @test norm(K_ff - K_ff2) / norm(K_ff2) <= 1.0e-5
    @test abs(v_f' * K_ff2 * v_f - (result)) / (result) <= 1.0e-5
    true
end
test()
nothing
end

module mparallelassembly_high_level_1
using FinEtools
using FinEtoolsMultithreading.Exports
using LinearAlgebra
using Test

function test()
    W = 1.1
    L = 12.0
    t = 0.32
    nl, nt, nw = 12, 33, 24
    nl, nt, nw = 2, 3, 2
    ntasks = 2

    fens, fes = H8block(L, W, t, nl, nw, nt)
    geom = NodalField(fens.xyz)
    psi = NodalField(fill(1.0, count(fens), 1))
    nl = collect(1:3)
    setebc!(psi, nl, true, ones(Int, length(nl)), 0.0)
    numberdofs!(psi)
    v_f = gathersysvec(psi)

    function createsubdomain(fessubset)
        FEMMBase(IntegDomain(fessubset, GaussRule(3, 2)))
    end

    function matrixcomputation!(femm, assembler)
        bilform_diffusion(femm, assembler, geom, psi, DataCache(Matrix(1.0 * LinearAlgebra.I(3))))
    end

    K = parallel_make_matrix(
        fes,
        psi,
        createsubdomain,
        matrixcomputation!;
        ntasks=Threads.nthreads(),
        kind=:CSC,
    )

    K_ff = matrix_blocked_ff(K, nfreedofs(psi))
    result = abs(v_f' * K_ff * v_f)

    ass = SysmatAssemblerFFBlock(nfreedofs(psi))
    K_ff2 = bilform_diffusion(FEMMBase(IntegDomain(fes, GaussRule(3, 2))), ass, geom, psi, DataCache(Matrix(1.0 * LinearAlgebra.I(3))))
    @test norm(K_ff - K_ff2) / norm(K_ff2) <= 1.0e-5
    @test abs(v_f' * K_ff2 * v_f - (result)) / (result) <= 1.0e-5
    true
end
test()
nothing
end

module mparallelassembly_high_level_2
using FinEtools
using FinEtoolsMultithreading.Exports
using FinEtoolsMultithreading: domain_decomposition, parallel_matrix_assembly!, SysmatAssemblerSparsePatt, associate_pattern
using LinearAlgebra
using Test

function test()
    W = 1.1
    L = 12.0
    t = 0.32
    nl, nt, nw = 12, 33, 24
    nl, nt, nw = 2, 3, 2
    ntasks = 2

    fens, fes = H8block(L, W, t, nl, nw, nt)
    geom = NodalField(fens.xyz)
    psi = NodalField(fill(1.0, count(fens), 1))
    nl = collect(1:3)
    setebc!(psi, nl, true, ones(Int, length(nl)), 0.0)
    numberdofs!(psi)
    v_f = gathersysvec(psi)

    function createsubdomain(fessubset)
        FEMMBase(IntegDomain(fessubset, GaussRule(3, 2)))
    end

    function matrixcomputation!(femm, assembler)
        bilform_diffusion(femm, assembler, geom, psi, DataCache(Matrix(1.0 * LinearAlgebra.I(3))))
    end

    n2e = FENodeToFEMap(fes.conn, count(fens))
    n2n = FENodeToNeighborsMap(n2e, fes.conn)
    K = sparse_symmetric_csc_pattern(psi.dofnums, nents(psi), n2n, zero(eltype(psi.values)))
    element_colors, unique_colors = element_coloring(fes, n2e)
    decomposition = domain_decomposition(fes, ntasks, element_colors, unique_colors, createsubdomain)
    assembler = SysmatAssemblerSparsePatt(0.0)
    associate_pattern(assembler, K)
    K = parallel_matrix_assembly!(
        assembler,
        decomposition,
        matrixcomputation!,
        ntasks
    )

    K_ff = matrix_blocked_ff(K, nfreedofs(psi))
    result = abs(v_f' * K_ff * v_f)

    ass = SysmatAssemblerFFBlock(nfreedofs(psi))
    K_ff2 = bilform_diffusion(FEMMBase(IntegDomain(fes, GaussRule(3, 2))), ass, geom, psi, DataCache(Matrix(1.0 * LinearAlgebra.I(3))))
    @test norm(K_ff - K_ff2) / norm(K_ff2) <= 1.0e-5
    @test abs(v_f' * K_ff2 * v_f - (result)) / (result) <= 1.0e-5
    true
end
test()
nothing
end

module mparallelassembly_high_level_3
using FinEtools
using FinEtoolsMultithreading.Exports
using FinEtoolsMultithreading: domain_decomposition, parallel_matrix_assembly!, SysmatAssemblerSparsePatt, associate_pattern
using LinearAlgebra
using Test

function test()
    W = 1.1
    L = 12.0
    t = 0.32
    nl, nt, nw = 12, 33, 24
    ntasks = 2

    fens, fes = H8block(L, W, t, nl, nw, nt)
    geom = NodalField(fens.xyz)
    psi = NodalField(fill(1.0, count(fens), 1))
    nl = collect(1:3)
    setebc!(psi, nl, true, ones(Int, length(nl)), 0.0)
    numberdofs!(psi)
    v_f = gathersysvec(psi)

    function createsubdomain(fessubset)
        FEMMBase(IntegDomain(fessubset, GaussRule(3, 2)))
    end

    function matrixcomputation!(femm, assembler)
        bilform_diffusion(femm, assembler, geom, psi, DataCache(Matrix(1.0 * LinearAlgebra.I(3))))
    end

    n2e = FENodeToFEMap(fes.conn, count(fens))
    n2n = FENodeToNeighborsMap(n2e, fes.conn)
    K = sparse_symmetric_csc_pattern(psi.dofnums, nents(psi), n2n, zero(eltype(psi.values)))
    element_colors, unique_colors = element_coloring(fes, n2e)
    decomposition = domain_decomposition(fes, ntasks, element_colors, unique_colors, createsubdomain)
    assembler = SysmatAssemblerSparsePatt(0.0)
    associate_pattern(assembler, K)
    K = parallel_matrix_assembly!(
        assembler,
        decomposition,
        matrixcomputation!,
        ntasks
    )

    K_ff = matrix_blocked_ff(K, nfreedofs(psi))
    result = abs(v_f' * K_ff * v_f)

    ass = SysmatAssemblerFFBlock(nfreedofs(psi))
    K_ff2 = bilform_diffusion(FEMMBase(IntegDomain(fes, GaussRule(3, 2))), ass, geom, psi, DataCache(Matrix(1.0 * LinearAlgebra.I(3))))
    @test norm(K_ff - K_ff2) / norm(K_ff2) <= 1.0e-5
    @test abs(v_f' * K_ff2 * v_f - (result)) / (result) <= 1.0e-5
    true
end
test()
nothing
end
