"""
    parallel_make_matrix(fes, u, crsubdom, matrixcomputation!, ntasks, kind = :CSC)

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
    crsubdom,
    matrixcomputation!;
    ntasks = Threads.nthreads(),
    kind = :CSC,
)
    n2e = FENodeToFEMap(fes.conn, nnodes(u))
    parallel_make_matrix(
        fes,
        u,
        n2e,
        crsubdom,
        matrixcomputation!,
        ntasks,
        kind
    )
end

"""
    parallel_make_matrix(
        fes,
        u,
        n2e,
        crsubdom,
        matrixcomputation!,
        ntasks,
        kind = :CSC,
    )

Assemble a sparse matrix.
"""
function parallel_make_matrix(
    fes,
    u,
    n2e,
    crsubdom,
    matrixcomputation!,
    ntasks,
    kind
)
    @assert kind in [:CSC, :CSR]
    n2n = FENodeToNeighborsMap(n2e, fes.conn)
    K = sparse_symmetric_csc_pattern(u.dofnums, nents(u), n2n, zero(eltype(u.values)))
    element_colors, unique_colors = element_coloring(fes, n2e)
    decomposition = domain_decomposition(fes, ntasks, element_colors, unique_colors, crsubdom)
    assembler = SysmatAssemblerSparsePatt(0.0)
    associate_pattern(assembler, K)
    parallel_matrix_assembly!(
        assembler,
        decomposition,
        matrixcomputation!,
        ntasks
    )
    return assembler._pattern
end
