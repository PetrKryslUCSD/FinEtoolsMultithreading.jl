"""
    parallel_make_matrix(
        fes,
        u,
        createsubd,
        matrixupdt!;
        ntasks = Threads.nthreads(),
        kind = :CSC,
    )

Assemble a sparse matrix.

Either a `:CSC` matrix or a `:CSR` matrix is created. We shall refer to this
matrix as a CSX matrix. The process is:

1. Construct the incidence relation node-to-neighbors.
2. Make the sparse pattern and create a sparse CSX matrix with all values zero.
3. Construct the incidence relation element-to-neighbors.
4. Compute element coloring.
5. Set up domain decomposition.
6. Compute and assemble the matrix entries.
"""
function parallel_make_matrix(
    fes,
    u,
    createsubd,
    matrixupdt!;
    ntasks = Threads.nthreads(),
    kind = :CSC,
)
    n2e = FENodeToFEMap(fes, nnodes(u))
    parallel_make_matrix(
        fes,
        u.dofnums,
        nalldofs(u),
        eltype(u.values),
        n2e,
        createsubd,
        matrixupdt!,
        ntasks,
        kind,
    )
end

"""
    parallel_make_matrix(
        fes,
        dofnums,
        ndofs,
        FT,
        n2e,
        createsubd,
        matrixupdt!,
        ntasks,
        kind,
    )

Assemble a sparse matrix.

1. Construct the incidence relation node-to-neighbors.
2. Make the sparse pattern and create a sparse CSX matrix with all values zero.
3. Construct the incidence relation element-to-neighbors.
4. Compute element coloring.
5. Set up domain decomposition.
6. Compute and assemble the matrix entries.
"""
function parallel_make_matrix(
    fes,
    dofnums,
    ndofs,
    FT,
    n2e,
    createsubd,
    matrixupdt!,
    ntasks,
    kind,
)
    @assert kind in [:CSC, :CSR]
    n2n = FENodeToNeighborsMap(n2e, fes)
    K_pattern = sparse_symmetric_csc_pattern(dofnums, ndofs, n2n, zero(FT))
    e2e = FEElemToNeighborsMap(n2e, fes)
    coloring = element_coloring(fes, e2e)
    decomposition = domain_decomposition(fes, coloring, createsubd, ntasks)
    return parallel_matrix_assembly!(
        SysmatAssemblerSparsePatt(K_pattern),
        decomposition,
        matrixupdt!,
    )
end

