module mscan001
using Test
using ThreadedScans
using FinEtoolsMultithreading: _scan!, psp_scan!
function test()
    a = [1, 3, 4, 2, 6, 3, 2, 1, 3, 4, 7, 8, 3]
    for i in 1:12
        a = vcat(a, a)
    end
    @show length(a)
    
    @info "ThreadedScans.scan!"
    b = deepcopy(a)
    ThreadedScans.scan!(+, b[1:5])
    ThreadedScans.scan!(+, b)
    
    @info "_scan!"
    c = deepcopy(a)
    _scan!(c)
    
    @test maximum(abs.(b - c)) == 0
    
    @info "psp_scan!"
    d = deepcopy(a)
    psp_scan!(d[1:5])
    d = deepcopy(a)
    psp_scan!(d)
    
    @test maximum(abs.(b - d)) == 0
    true
end
test()
nothing
end

module mparallelfenodetoelemmap
using FinEtools
using FinEtoolsMultithreading.Exports
using FinEtoolsMultithreading: SysmatAssemblerSparsePatt, SysmatAssemblerSparsePattwLookup, startassembly!, assemble!, makematrix!
using FinEtoolsMultithreading.FENodeToFEMapModule: FENodeToFEMapThr
using LinearAlgebra
using Test

function test()
    W = 1.1
    L = 12.0
    t = 0.32
    nl, nt, nw = 15, 13, 12
    ntasks = 2

    fens, fes = H8block(L, W, t, nl, nw, nt)

    n2e = FinEtools.FENodeToFEMap(fes.conn, count(fens))
    # n2e = FinEtools.FENodeToFEMap(fes.conn, count(fens))
    # n2ethr = FENodeToFEMapThr(fes, count(fens))
    n2ethr = FENodeToFEMapThr(fes, count(fens))
    for i in eachindex(n2e.map)
        @test n2e.map[i] == n2ethr.map[i]
    end
    # @info "On $(Threads.nthreads()) threads"
    true
end
test()
nothing
end

module mparallelassembly_assembler_w_lup_2
using FinEtools
using FinEtoolsMultithreading.Exports
using FinEtoolsMultithreading: SysmatAssemblerSparsePatt, SysmatAssemblerSparsePattwLookup, startassembly!, assemble!, makematrix!
using FinEtoolsMultithreading.FENodeToFEMapModule: FENodeToFEMapThr
using LinearAlgebra
using Test

function test()
    W = 1.1
    L = 12.0
    t = 0.32
    nl, nt, nw = 12, 33, 24
    nl, nt, nw = 5, 3, 4
    ntasks = 2

    fens, fes = H8block(L, W, t, nl, nw, nt)
    geom = NodalField(fens.xyz)
    psi = NodalField(fill(1.0, count(fens), 1))
    nl = collect(1:3)
    setebc!(psi, nl, true, ones(Int, length(nl)), 0.0)
    numberdofs!(psi)
    v_f = gathersysvec(psi)

    femm = FEMMBase(IntegDomain(fes, GaussRule(3, 2)))
    n2e = FENodeToFEMapThr(fes, count(fens))
    n2n = FENodeToNeighborsMap(n2e, fes.conn)
    K = csc_symmetric_pattern(psi.dofnums, nalldofs(psi), n2n, eltype(psi.values))
    K = csc_symmetric_pattern(psi.dofnums, nalldofs(psi), n2n)
    
    c = [i for i in fes.conn[end]]
    z = zeros(8, 8); r = psi.dofnums[c, 1]
    
    assembler =  SysmatAssemblerSparsePattwLookup(K)
    startassembly!(assembler, 8, 8, 1000, nalldofs(psi), nalldofs(psi))
    assemble!(assembler, z, r, r)
    assembler =  SysmatAssemblerSparsePatt(K)
    startassembly!(assembler, 8, 8, 1000, nalldofs(psi), nalldofs(psi))
    assemble!(assembler, z, r, r)
#    @code_warntype  assemble!(assembler, zeros(8, 8), psi.dofnums[c, 1], psi.dofnums[c, 1])

    true
end
test()
nothing
end

module mparallelassembly_assembler_w_lup_1
using FinEtools
using FinEtoolsMultithreading.Exports
using FinEtoolsMultithreading: SysmatAssemblerSparsePattwLookup, startassembly!, assemble!, makematrix!
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
    K = csc_symmetric_pattern(psi.dofnums, nalldofs(psi), n2n, eltype(psi.values))
    assmblr = SysmatAssemblerSparsePattwLookup(K)
    K = bilform_diffusion(femm, assmblr, geom, psi, DataCache(Matrix(1.0 * LinearAlgebra.I(3))))
    
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

module mmmmaps1
using FinEtools
using FinEtoolsMultithreading.Exports
using Test
function test(n = 2)
    W, L, H = 3.5, 7.1, 9.3
    fens, fes = H8block(W, L, H, n, n, n)
 
    n2e = FENodeToFEMap(fes.conn, count(fens))
    e2e = FElemToNeighborsMap(n2e, fes.conn)
    
    found = true
    for i in eachindex(n2e.map)
        for k in n2e.map[i]
            for m in n2e.map[i]
                (k != m) && ( # self-references are excluded
                found = found && ((k in e2e.map[m]) && (m in e2e.map[k]))
                )
            end
        end
    end
    @test found

    e2e = FElemToNeighborsMap(n2e, fes)
    
    found = true
    for i in eachindex(n2e.map)
        for k in n2e.map[i]
            for m in n2e.map[i]
                (k != m) && ( # self-references are excluded
                found = found && ((k in e2e.map[m]) && (m in e2e.map[k]))
                )
            end
        end
    end
    @test found
end
test(2)
test(17)
test(18)
nothing
end

module mmmmaps2
using FinEtools
using FinEtoolsMultithreading.Exports
using Test
function test(n = 20)
    W, L, H = 3.5, 7.1, 9.3
    fens, fes = H8block(W, L, H, n, n, n)
   
    # println("nalldofs(u) = $(nalldofs(u))").#
    n2e = FENodeToFEMap(fes.conn, count(fens))
    n2n = FENodeToNeighborsMap(n2e, fes.conn)
    
    found = true
    for i in eachindex(fes.conn)
        for k in fes.conn[i]
            for m in fes.conn[i]
                (k != m) && ( # exclude self reference
                found = found && ((k in n2n.map[m]) && (m in n2n.map[k]))
                )
            end
        end
    end
    @test found

    n2n = FENodeToNeighborsMap(n2e, fes)
    
    found = true
    for i in eachindex(fes.conn)
        for k in fes.conn[i]
            for m in fes.conn[i]
                (k != m) && ( # exclude self reference
                found = found && ((k in n2n.map[m]) && (m in n2n.map[k]))
                )
            end
        end
    end
    @test found
end
test(13)
test(17)
test(18)
nothing
end

module mparallelassembly_assembler_1
using FinEtools
using FinEtoolsMultithreading.Exports
using FinEtoolsMultithreading: SysmatAssemblerSparsePatt, startassembly!, assemble!, makematrix!
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
    K = csc_symmetric_pattern(psi.dofnums, nalldofs(psi), n2n, eltype(psi.values))
    assmblr = SysmatAssemblerSparsePatt(K)
    K = bilform_diffusion(femm, assmblr, geom, psi, DataCache(Matrix(1.0 * LinearAlgebra.I(3))))
    
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
using FinEtoolsMultithreading: decompose, parallel_matrix_assembly!, SysmatAssemblerSparsePatt
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
    K_pattern = csc_symmetric_pattern(psi.dofnums, nalldofs(psi), n2n, eltype(psi.values))
    coloring = element_coloring(fes, n2e)
    decomposition = decompose(fes, coloring, createsubdomain, ntasks)
    K = parallel_matrix_assembly!(
        SysmatAssemblerSparsePatt(K_pattern),
        decomposition,
        matrixcomputation!
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
using FinEtoolsMultithreading: decompose, parallel_matrix_assembly!, SysmatAssemblerSparsePatt
using LinearAlgebra
using Test
function test_coloring(coloring, n2e)
    element_colors, unique_colors = coloring
    @assert norm(sort(unique(element_colors)) - sort(unique_colors)) == 0
    for k in eachindex(n2e.map)
        nc = element_colors[n2e.map[k]]
        @assert length(nc) == length(unique(nc))
        @assert norm(sort(nc) - sort(unique(nc))) == 0
    end
end
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
    K_pattern = csc_symmetric_pattern(psi.dofnums, nalldofs(psi), n2n, eltype(psi.values))
    coloring = element_coloring(fes, n2e)
    test_coloring(coloring, n2e)

    decomposition = decompose(fes, coloring, createsubdomain, ntasks)
    K = parallel_matrix_assembly!(
        SysmatAssemblerSparsePatt(K_pattern),
        decomposition,
        matrixcomputation!
    )

    K_ff = matrix_blocked_ff(K, nfreedofs(psi))
    result = abs(v_f' * K_ff * v_f)

    ass = SysmatAssemblerFFBlock(nfreedofs(psi))
    K_ff2 = bilform_diffusion(FEMMBase(IntegDomain(fes, GaussRule(3, 2))), ass, geom, psi, DataCache(Matrix(1.0 * LinearAlgebra.I(3))))
    if !(norm(K_ff - K_ff2) / norm(K_ff2) <= 1.0e-5)
        @show K_ff
        @show K_ff2
    end
    @test norm(K_ff - K_ff2) / norm(K_ff2) <= 1.0e-5
    @test abs(v_f' * K_ff2 * v_f - (result)) / (result) <= 1.0e-5
    true
end
test()
nothing
end

module mparallelassembly_high_level_4
using FinEtools
using FinEtoolsMultithreading.Exports
using FinEtoolsMultithreading: decompose, parallel_matrix_assembly!, SysmatAssemblerSparsePatt
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
    K = parallel_make_matrix(
        fes,
        psi,
        createsubdomain,
        matrixcomputation!;
        ntasks,
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

module mparallelassembly_high_level_5
using FinEtools
using FinEtoolsMultithreading.Exports
using FinEtoolsMultithreading: decompose, parallel_matrix_assembly!, SysmatAssemblerSparsePatt
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
    K = parallel_make_matrix(
        fes,
        psi,
        createsubdomain,
        matrixcomputation!;
        ntasks,
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

module mparallelassembly_high_level_6
using FinEtools
using FinEtoolsMultithreading.Exports
using FinEtoolsMultithreading: decompose, parallel_matrix_assembly!, SysmatAssemblerSparsePatt
using LinearAlgebra
using Test

function test()
    W = 1.1
    L = 12.0
    t = 0.32
    nl, nt, nw = 12, 33, 24
    ntasks = 1

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
        matrixcomputation!
    )


    # n2e = FENodeToFEMap(fes, nnodes(psi))
    # n2n = FENodeToNeighborsMap(n2e, fes)
    # K_pattern = csc_symmetric_pattern(psi.dofnums, nalldofs(psi), n2n, zero(Float64))
    # coloring = element_coloring(fes, n2e)
    # decomposition = decompose(fes, coloring, createsubdomain, ntasks)
    # K = parallel_matrix_assembly!(
    #     SysmatAssemblerSparsePatt(K_pattern),
    #     decomposition,
    #     matrixcomputation!
    # )

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
