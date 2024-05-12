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

1. Construct the incidence relation node-to-elements.
2. Construct the incidence relation node-to-neighbors.
3. Make the sparse pattern and create a sparse CSX matrix with all values zero.
4. Construct the incidence relation element-to-neighbors.
5. Compute element coloring.
6. Set up domain decomposition.
7. Compute and assemble the matrix entries.

Here, `createsubd` could be
```
(fessubset) -> FEMMAcoust(IntegDomain(fessubset, GaussRule(3, 2)), material)
```
and `matrixupdt!` could be
```
(femm, assmblr) -> acousticstiffness(femm, assmblr, geom, P)
```

"""
function parallel_make_matrix(
    fes,
    u,
    createsubd,
    matrixupdt!;
    ntasks = Threads.nthreads(),
    kind = :CSC,
)
    n2e = FENodeToFEMapThr(fes, nnodes(u))
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
    fes, dofnums, ndofs, FT, n2e, # data
    createsubd, matrixupdt!,      # functions
    ntasks,                       # how many threads?
    kind
)
    @assert kind in [:CSC, :CSR]
    n2n = FENodeToNeighborsMap(n2e, fes) # ALG 1
    matrix = csc_symmetric_pattern(dofnums, # ALG 2
                                 ndofs, n2n, FT)
    e2e = FElemToNeighborsMap(n2e, fes) # ALG 3
    coloring = element_coloring(fes, e2e, ntasks) # ALG 4
    decomposition = decompose(fes, coloring, # ALG 5
                              createsubd, ntasks)
    return parallel_matrix_assembly!(  # ALG 6
        SysmatAssemblerSparsePatt(matrix),
        decomposition,
        matrixupdt!,
    )
end

