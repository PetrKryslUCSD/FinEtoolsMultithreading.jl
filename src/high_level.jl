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
    return parallel_make_matrix(
        fes,
        u,
        n2e,
        crsubdom,
        matrixcomputation!,
        ntasks,
        kind,
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
    assblr = fill_assembler(fes, u, crsubdom, matrixcomputation!, ntasks)
    n2n = FENodeToNeighborsMap(n2e, fes.conn)
    K = sparse_symmetric_zero(u, n2n, kind)
    return add_to_matrix!(K, assblr)
end
