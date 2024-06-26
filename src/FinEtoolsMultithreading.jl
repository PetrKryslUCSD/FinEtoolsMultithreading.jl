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
using ECLGraphColor

include("utilities.jl")
include("prefixsum.jl")
include("FENodeToNeighborsMapModule.jl")
using .FENodeToNeighborsMapModule: FENodeToNeighborsMap
include("FElemToNeighborsMapModule.jl")
using .FElemToNeighborsMapModule: FElemToNeighborsMap
include("FENodeToFEMapModule.jl")
include("element_coloring.jl")
include("parallel_element_coloring.jl")
include("sparsity_pattern.jl")
include("parallel_assembly.jl")
include("domain_decomposition.jl")
using .FENodeToFEMapModule: FENodeToFEMapThr
include("high_level.jl")


module Exports
# The high level interface
using ..FinEtoolsMultithreading: parallel_make_matrix
export parallel_make_matrix
# These three routines give access to intermediate steps
using ..FinEtoolsMultithreading: parallel_matrix_assembly!
export parallel_matrix_assembly!
using ..FinEtoolsMultithreading: csc_symmetric_pattern
export csc_symmetric_pattern
using ..FENodeToNeighborsMapModule: FENodeToNeighborsMap
# Exported: type for maps from nodes to connected nodes
export FENodeToNeighborsMap
using ..FElemToNeighborsMapModule: FElemToNeighborsMap
# Exported: type for maps from nodes to connected nodes
export FElemToNeighborsMap
using ..FENodeToFEMapModule: FENodeToFEMapThr
# Exported: type for maps from nodes to elements computed on multiple threads
export FENodeToFEMapThr
export element_coloring
end

# Enable LSP look up in test modules
if false
    include("../test/runtests.jl")
end

end # module FinEtoolsMultithreading
