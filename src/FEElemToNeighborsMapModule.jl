"""
    FEElemToNeighborsMapModule  

Module to construct a map from finite elements to elements that share nodes with
them (i.e. that are connected to them through the nodes).
"""
module FEElemToNeighborsMapModule

__precompile__(true)

using FinEtools.FESetModule: AbstractFESet
using FinEtools.FENodeToFEMapModule: FENodeToFEMap

function _unique_elem_neighbors(n2emap, conn1)
    neighbors = fill(zero(eltype(conn1)), 0)
    @inbounds for k in conn1
        for j in n2emap[k]
            push!(neighbors, j)
        end
    end
    return unique!(sort!(neighbors))
end

function _e2e_map(n2e, conn)
    T = typeof(n2e.map[1])
    map = Array{T}(undef, length(conn))
    Base.Threads.@threads for i in eachindex(map) # run this in PARALLEL
        map[i] = _unique_elem_neighbors(n2e.map, conn[i])
    end
    return map
end

"""
    FEElemToNeighborsMap

Map from finite elements to the elements that are connected to them by finite
nodes.

!!! note

    Self references are included (an element is connected to itself).
"""
struct FEElemToNeighborsMap{IT}
    # Map as a vector of vectors.
    map::Vector{Vector{IT}}
end

"""
    FEElemToNeighborsMap(
        n2e::N2EMAP,
        conn::Vector{NTuple{N,IT}},
    ) where {N2EMAP<:FENodeToFEMap,N,IT<:Integer}

Map from finite elements to the elements that are connected by finite nodes.

- `conns` = connectivities as a vector of tuples
- `nmax` = largest possible node number
"""
function FEElemToNeighborsMap(
    n2e::N2EMAP,
    conn::Vector{NTuple{N,IT}},
) where {N2EMAP<:FENodeToFEMap,N,IT<:Integer}
    return FEElemToNeighborsMap(_e2e_map(n2e, conn))
end

"""
    FEElemToNeighborsMap(
        n2e::N2EMAP,
        fes::FE,
    ) where {N2EMAP<:FENodeToFEMap,FE<:AbstractFESet}

Map from finite elements to the elements that are connected by finite nodes.

Convenience constructor.
"""
function FEElemToNeighborsMap(
    n2e::N2EMAP,
    fes::FE,
) where {N2EMAP<:FENodeToFEMap,FE<:AbstractFESet}
    return FEElemToNeighborsMap(_e2e_map(n2e, fes.conn))
end

end
