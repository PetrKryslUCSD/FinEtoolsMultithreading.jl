"""
    FinEtoolsMultithreading

Package for parallel finite element computing.
"""
module FinEtoolsMultithreading

using SparseArrays
using ChunkSplitters
using SparseMatricesCSR
using FinEtools
using ThreadedScans


include("utilities.jl")
include("parallel_sorting.jl")
using .PQuickSort
include("sparsity_pattern.jl")
include("parallel_assembly.jl")
include("add_to_sparse.jl")
include("domain_decomposition.jl")
include("high_level.jl")

module Exports
# The high level interface
using ..FinEtoolsMultithreading: parallel_make_matrix
export parallel_make_matrix
# These three routines give access to intermediate steps
using ..FinEtoolsMultithreading: fill_assembler
export fill_assembler
using ..FinEtoolsMultithreading: sparse_symmetric_zero
export sparse_symmetric_zero
using ..FinEtoolsMultithreading: add_to_matrix!
export add_to_matrix!
end

# Enable LSP look up in test modules
if false
    include("../test/runtests.jl")
end

end # module FinEtoolsMultithreading
