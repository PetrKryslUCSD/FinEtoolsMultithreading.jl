"""
    FElemToNeighborsMapModule  

Module to construct a map from finite elements to elements that share nodes with
them (i.e. that are connected to them through the nodes).
"""
module FElemToNeighborsMapModule

__precompile__(true)

using FinEtools.FESetModule: AbstractFESet
using FinEtools.FENodeToFEMapModule: FENodeToFEMap

function _unique_elem_neighbors(self, n2emap, conn1)
    len = 0
    @inbounds for k in conn1
        len += (length(n2emap[k]) - 1) # we are not adding self-references
    end
    neighbors = fill(zero(eltype(conn1)), len)
    p = 1
    @inbounds for k in conn1
        for j in n2emap[k]
            if j != self # we are not adding self-reference
                neighbors[p] = j; p += 1
            end
        end
    end
    return unique!(sort!(neighbors))
end

function _e2e_map(n2e, conn)
    T = typeof(n2e.map[1])
    map = Array{T}(undef, length(conn))
    Base.Threads.@threads for i in eachindex(map) # run this in PARALLEL
        map[i] = _unique_elem_neighbors(i, n2e.map, conn[i])
    end
    return map
end

"""
    FElemToNeighborsMap

Map from finite elements to the elements that are connected to them by finite
nodes.

!!! note

    Self references are included (an element is connected to itself).
"""
struct FElemToNeighborsMap{IT}
    # Map as a vector of vectors.
    map::Vector{Vector{IT}}
end

"""
    FElemToNeighborsMap(
        n2e::N2EMAP,
        conn::Vector{NTuple{N,IT}},
    ) where {N2EMAP<:FENodeToFEMap,N,IT<:Integer}

Map from finite elements to the elements that are connected by finite nodes.

- `conns` = connectivities as a vector of tuples
- `nmax` = largest possible node number
"""
function FElemToNeighborsMap(
    n2e::N2EMAP,
    conn::Vector{NTuple{N,IT}},
) where {N2EMAP<:FENodeToFEMap,N,IT<:Integer}
    return FElemToNeighborsMap(_e2e_map(n2e, conn))
end

"""
    FElemToNeighborsMap(
        n2e::N2EMAP,
        fes::FE,
    ) where {N2EMAP<:FENodeToFEMap,FE<:AbstractFESet}

Map from finite elements to the elements that are connected by finite nodes.

Convenience constructor.
"""
function FElemToNeighborsMap(
    n2e::N2EMAP,
    fes::FE,
) where {N2EMAP<:FENodeToFEMap,FE<:AbstractFESet}
    return FElemToNeighborsMap(_e2e_map(n2e, fes.conn))
end

end
