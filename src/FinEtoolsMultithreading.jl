"""
    FinEtoolsMultithreading

Package for parallel finite element computing.
"""
module FinEtoolsMultithreading

using SparseArrays
using ChunkSplitters
using SparseMatricesCSR
using LinearAlgebra
using FinEtools
using ThreadedScans

include("utilities.jl")
include("FENodeToNeighborsMapModule.jl")
using .FENodeToNeighborsMapModule: FENodeToNeighborsMap
include("FEElemToNeighborsMapModule.jl")
using .FEElemToNeighborsMapModule: FEElemToNeighborsMap
include("element_coloring.jl")
include("sparsity_pattern.jl")
include("parallel_assembly.jl")
include("domain_decomposition.jl")
include("high_level.jl")


module Exports
# The high level interface
using ..FinEtoolsMultithreading: parallel_make_matrix
export parallel_make_matrix
# These three routines give access to intermediate steps
using ..FinEtoolsMultithreading: parallel_matrix_assembly!
export parallel_matrix_assembly!
using ..FinEtoolsMultithreading: sparse_symmetric_csc_pattern
export sparse_symmetric_csc_pattern
using ..FENodeToNeighborsMapModule: FENodeToNeighborsMap
# Exported: type for maps from nodes to connected nodes
export FENodeToNeighborsMap
using ..FEElemToNeighborsMapModule: FEElemToNeighborsMap
# Exported: type for maps from nodes to connected nodes
export FEElemToNeighborsMap
end

# Enable LSP look up in test modules
if false
    include("../test/runtests.jl")
end

end # module FinEtoolsMultithreading
