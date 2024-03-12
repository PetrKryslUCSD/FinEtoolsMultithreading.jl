
function parallel_make_csc_matrix(fes, u, crsubdom, matrixcomputation!, ntasks)
    femms = subdomainfemms(fes, ntasks, crsubdom)
    assembler = make_assembler(femms, SysmatAssemblerSparse, u)
    start_assembler!(assembler)
    assemblers = make_task_assemblers(femms, assembler, SysmatAssemblerSparse, u)
    parallel_matrix_assembly(femms, assemblers, matrixcomputation!)
    K = csc_matrix_pattern(fes, u)
    return add_to_matrix!(K, assembler)
end