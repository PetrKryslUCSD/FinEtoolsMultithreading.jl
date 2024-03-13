"""
    ParFEM

Package for parallel finite element computing.
"""
module ParFEM

using SparseArrays
using ChunkSplitters
using SparseMatricesCSR
using FinEtools

include("sparsity_pattern.jl")
include("add_to_sparse.jl")
include("parallel_assembly.jl")
include("domain_decomposition.jl")
include("high_level.jl")

# Enable LSP look up in test modules
if false
    include("../test/runtests.jl")
end

end # module ParFEM
