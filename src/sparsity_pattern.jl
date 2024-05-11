function _store_dofs!(rowval, n, nghbrs, dofnums, colptr)
    s1 = colptr[dofnums[n, 1]]
    p = 0
    k = n # The node itself needs to be considered
    @inbounds for d in axes(dofnums, 2)
        rowval[s1+p] = dofnums[k, d]; p += 1
    end
    @inbounds for k in nghbrs
        for d in axes(dofnums, 2)
            rowval[s1+p] = dofnums[k, d]; p += 1
        end
    end
    bl = p # block length
    sort!(@view(rowval[s1:s1+bl-1]))
    @inbounds for d in 2:size(dofnums, 2)
        s = colptr[dofnums[n, d]]
        copy!(@view(rowval[s:s+bl-1]), @view(rowval[s1:s1+bl-1]))
    end
    return nothing
end

function _row_blocks(IT, map, dofnums)
    nd = size(dofnums, 2)
    total_dofs = length(map) * nd
    lengths = Vector{IT}(undef, total_dofs + 1)
    lengths[1] = 1
    @inbounds Threads.@threads for k in eachindex(map)
        kl = (length(map[k]) + 1) * nd # +1 for the node itself
        for d in axes(dofnums, 2)
            j = dofnums[k, d]
            lengths[j+1] = kl
        end
    end
    return lengths
end

function _acc_start_ptr!(s)
    len = length(s)
    @inbounds for k = 1:len-1
        s[k+1] += s[k]
    end
    s
end

function _csc_arrays(IT, FT, map, dofnums)
    # First we create an array of the lengths of the dof blocks
    colptr = _row_blocks(IT, map, dofnums)
    # Now we start overwriting the "lengths" array with the starts
    ThreadedScans.scan!(+, colptr) # equivalent to _acc_start_ptr!(start)
    sumlen = colptr[end] - 1
    rowval = Vector{IT}(undef, sumlen) # This will get filled in later
    nzval = _zeros_via_calloc(FT, sumlen) # This needs to be initialized for future accumulation
    return colptr, rowval, nzval
end

"""
    csc_symmetric_pattern(
        dofnums::Array{IT,2},
        nrowscols,
        n2n,
        FT = Float64,
    ) where {IT<:Integer}
    
Create symmetric sparse zero matrix (sparsity pattern).

Uses the following data structures:
```
    n2n = FENodeToNeighborsMap(n2e, fes.conn)
```
"""
function csc_symmetric_pattern(
    dofnums::Array{IT,2},
    nrowscols,
    n2n,
    FT = Float64,
) where {IT<:Integer}
    @assert length(n2n.map) == size(dofnums, 1)
    colptr, rowval, nzval = _csc_arrays(IT, FT, n2n.map, dofnums)
    Base.Threads.@threads for n in axes(dofnums, 1)
        _store_dofs!(rowval, n, n2n.map[n], dofnums, colptr)
    end 
    return SparseMatrixCSC(nrowscols, nrowscols, colptr, rowval, nzval)
end
