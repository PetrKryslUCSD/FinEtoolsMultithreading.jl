"""
    ParFEM

Package for parallel finite element computing.
"""
module ParFEM

using SparseArrays
using ChunkSplitters
using FinEtools

include("sparsity_pattern.jl")
include("add_to_sparse.jl")
include("parallel_assembly.jl")
include("domain_decomposition.jl")
include("high_level.jl")

end # module ParFEM
