"""
    parallel_make_matrix(fes, u, crsubdom, matrixcomputation!, ntasks, kind = :CSC)

Assemble a sparse matrix.
"""
function fill_assembler(fes, u, crsubdom, matrixcomputation!, ntasks)
    femms = subdomainfemms(fes, ntasks, crsubdom)
    assembler = make_assembler(femms, SysmatAssemblerSparse, u)
    start_assembler!(assembler)
    assemblers = make_task_assemblers(femms, assembler, SysmatAssemblerSparse, u)
    parallel_matrix_assembly(femms, assemblers, matrixcomputation!)
    return assembler
end

"""
    make_pattern_and_matrix(fes, u, kind = :CSC)

Make a sparse matrix from a sparsity pattern.
"""
function make_pattern_and_matrix(fes, u, kind = :CSC)
    @assert kind in [:CSC, :CSR]
    if kind == :CSC
        K = csc_matrix(fes, u)
    elseif kind == :CSR
        K = csr_matrix(fes, u)
    end
    return K
end

"""
    parallel_make_matrix(fes, u, crsubdom, matrixcomputation!, ntasks, kind = :CSC)

Assemble a sparse matrix.

Either a `:CSC` matrix or a `:CSR` matrix are created.

1. Compute the matrix entries as a COO sparse format.
2. Make the sparse pattern and create a sparse CSX matrix with all values zero.
3. Add the values from the COO list to the CSX sparse matrix.
"""
function parallel_make_matrix(fes, u, crsubdom, matrixcomputation!, ntasks, kind = :CSC)
    @assert kind in [:CSC, :CSR]
    assembler = fill_assembler(fes, u, crsubdom, matrixcomputation!, ntasks)
    K = make_pattern_and_matrix(fes, u, kind)
    return add_to_matrix!(K, assembler)
end
