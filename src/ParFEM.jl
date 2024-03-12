"""
    ParFEM

Package for parallel finite element computing.
"""
module ParFEM

using SparseArrays
using FinEtools

include("sparsity_pattern.jl")
include("add_to_sparse.jl")
include("parallel_matrix_assembly.jl")

end # module ParFEM
