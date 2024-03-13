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
using ..ParFEM: make_pattern
export fill_assembler, make_pattern
using ..ParFEM: add_to_matrix!
export add_to_matrix!
end

# Enable LSP look up in test modules
if false
    include("../test/runtests.jl")
end

end # module ParFEM
