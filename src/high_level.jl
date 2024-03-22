"""
    parallel_make_matrix(
        fes,
        u,
        createsubdomain,
        matrixcomputation!;
        ntasks = Threads.nthreads(),
        kind = :CSC,
    )

Assemble a sparse matrix.

Either a `:CSC` matrix or a `:CSR` matrix is created. We shall refer to this
matrix as a CSX matrix. The process is:

1. Compute the matrix entries as a COO sparse format;
2. Construct some incidence relation, and use them to
3. Make the sparse pattern and create a sparse CSX matrix with all values zero;
4. Add the values from the COO list to the CSX sparse matrix.
"""
function parallel_make_matrix(
    fes,
    u,
    createsubdomain,
    matrixcomputation!;
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
        createsubdomain, 
        matrixcomputation!, 
        ntasks, 
        kind)
end

"""
    parallel_make_matrix(
        fes,
        u,
        n2e,
        createsubdomain,
        matrixcomputation!,
        ntasks,
        kind
    )

Assemble a sparse matrix.

1. Make the sparse pattern and create a sparse CSX matrix with all values zero;
1. Compute the matrix entries as a COO sparse format;
2. Construct some incidence relation, and use them to
4. Add the values from the COO list to the CSX sparse matrix.
"""
function parallel_make_matrix(
    fes,
    dofnums,
    ndofs,
    FT,
    n2e,
    createsubdomain,
    matrixcomputation!,
    ntasks,
    kind,
)
    @assert kind in [:CSC, :CSR]
    n2n = FENodeToNeighborsMap(n2e, fes)
    K_pattern = sparse_symmetric_csc_pattern(dofnums, ndofs, n2n, zero(FT))
    coloring = element_coloring(fes, n2e)
    decomposition = domain_decomposition(fes, coloring, createsubdomain, ntasks)
    K = parallel_matrix_assembly!(SysmatAssemblerSparsePatt(K_pattern), decomposition, matrixcomputation!, ntasks)
    return K
end
