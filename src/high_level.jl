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
    make_pattern(fes, u, kind = :CSC)

Make sparsity pattern.
"""
function make_pattern(fes, u, kind = :CSC)
    @assert kind in [:CSC, :CSR]
    if kind == :CSC
        K = csc_matrix_pattern(fes, u)
    elseif kind == :CSR
        K = csr_matrix_pattern(fes, u)
    end
    return K
end

"""
    parallel_make_matrix(fes, u, crsubdom, matrixcomputation!, ntasks, kind = :CSC)

Assemble a sparse matrix.
"""
function parallel_make_matrix(fes, u, crsubdom, matrixcomputation!, ntasks, kind = :CSC)
    @assert kind in [:CSC, :CSR]
    assembler = fill_assembler(fes, u, crsubdom, matrixcomputation!, ntasks)
    if kind == :CSC
        K = csc_matrix_pattern(fes, u)
    elseif kind == :CSR
        K = csr_matrix_pattern(fes, u)
    end
    return add_to_matrix!(K, assembler)
end
