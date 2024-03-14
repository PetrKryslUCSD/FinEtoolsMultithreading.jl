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
include("parallel_assembly.jl")
include("add_to_sparse.jl")
include("domain_decomposition.jl")
include("high_level.jl")

module Exports
using ..ParFEM: parallel_make_matrix
export parallel_make_matrix
using ..ParFEM: fill_assembler
using ..ParFEM: make_pattern_and_matrix
export fill_assembler, make_pattern_and_matrix
using ..ParFEM: add_to_matrix!
export add_to_matrix!
using ..ParFEM: sparsity_pattern_symmetric
export sparsity_pattern_symmetric
end

# Enable LSP look up in test modules
if false
    include("../test/runtests.jl")
end

end # module ParFEM
