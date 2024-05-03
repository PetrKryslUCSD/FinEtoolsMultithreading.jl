"""
FENodeToFEMapModule  

Module to construct a map from finite element nodes to elements connected to them.
This map is constructed in parallel.1
"""
module FENodeToFEMapModule

__precompile__(true)

using FinEtools.FESetModule: AbstractFESet, subset, count
using FinEtools
using ChunkSplitters

function _merge_maps!(maps, rnge)
    for m in 2:length(maps)
        for i in rnge
            append!(maps[1][i], maps[m][i])
        end
    end
    return maps
end

function FENodeToFEMapThr(fes::FE, nmax::IT, ntasks = Threads.nthreads()) where {FE<:AbstractFESet,IT<:Integer}
    chks = chunks(1:count(fes), ntasks)
    maps = Vector{Vector{Vector{IT}}}(undef, length(chks))
    Base.Threads.@threads for c in chks
        maps[c[2]] = FinEtools.FENodeToFEMapModule._makemap(fes.conn, c[1], nmax)
    end       
    chks = chunks(1:nmax, ntasks)
    Base.Threads.@threads for c in chks
        _merge_maps!(maps, c[1])
    end    
    return FinEtools.FENodeToFEMapModule.FENodeToFEMap(maps[1])
end

end
