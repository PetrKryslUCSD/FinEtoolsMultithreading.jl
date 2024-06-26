"""
    FENodeToNeighborsMapModule  

Module to construct a map from finite element nodes to nodes connected to them.
"""
module FENodeToNeighborsMapModule

__precompile__(true)

using FinEtools.FESetModule: AbstractFESet
using FinEtools.FENodeToFEMapModule: FENodeToFEMap

function _nneighbors(self, ellist, conn, npe)
    totn = length(ellist) * (npe - 1) # we are not adding self-reference
    nodes = fill(zero(eltype(conn[1])), totn)
    p = 1
    @inbounds for i in ellist
        for k in conn[i]
            if k != self # we are not adding self-reference
                nodes[p] = k; p += 1
            end
        end
    end
    return unique!(sort!(nodes))
end

function _make_map(n2e, conn)
    npe = length(conn[1])
    map = Vector{eltype(n2e.map)}(undef, length(n2e.map))
    Base.Threads.@threads for i in eachindex(n2e.map) # run this in PARALLEL
        map[i] = _nneighbors(i, n2e.map[i], conn, npe)
    end
    return map
end

"""
    FENodeToNeighborsMap

Map from finite element nodes to the nodes that are connected to them by finite
elements.

!!! note

    Self references are excluded (a node is its own neighbor, 
    but it is not included in this data structure as such). 
"""
struct FENodeToNeighborsMap{IT}
    # Map as a vector of vectors.
    map::Vector{Vector{IT}}
end

"""
    FENodeToNeighborsMap(
        n2e::N2EMAP,
        conn::Vector{NTuple{N,IT}},
    ) where {N2EMAP<:FENodeToFEMap,N,IT<:Integer}

Map from finite element nodes to the nodes connected to them by elements.

- `conns` = connectivities as a vector of tuples
- `nmax` = largest possible node number
"""
function FENodeToNeighborsMap(
    n2e::N2EMAP,
    conn::Vector{NTuple{N,IT}},
) where {N2EMAP<:FENodeToFEMap,N,IT<:Integer}
    return FENodeToNeighborsMap(_make_map(n2e, conn))
end

"""
    FENodeToNeighborsMap(
        n2e::N2EMAP,
        fes::FE,
    ) where {N2EMAP<:FENodeToFEMap,FE<:AbstractFESet}

Map from finite element nodes to the nodes connected to them by elements.

Convenience constructor.
"""
function FENodeToNeighborsMap(
    n2e::N2EMAP,
    fes::FE,
) where {N2EMAP<:FENodeToFEMap,FE<:AbstractFESet}
    return FENodeToNeighborsMap(_make_map(n2e, fes.conn))
end

end
