"""
    parallel_make_matrix(fes, u, crsubdom, matrixcomputation!, ntasks, kind = :CSC)

Assemble a sparse matrix.
"""
function parallel_make_matrix(fes, u, crsubdom, matrixcomputation!, ntasks, kind = :CSC)
    @assert kind in [:CSC, :CSR]
    femms = subdomainfemms(fes, ntasks, crsubdom)
    assembler = make_assembler(femms, SysmatAssemblerSparse, u)
    start_assembler!(assembler)
    assemblers = make_task_assemblers(femms, assembler, SysmatAssemblerSparse, u)
    parallel_matrix_assembly(femms, assemblers, matrixcomputation!)
    if kind == :CSC
        K = csc_matrix_pattern(fes, u)
    elseif kind == :CSR
        K = csr_matrix_pattern(fes, u)
    end
    return add_to_matrix!(K, assembler)
end
