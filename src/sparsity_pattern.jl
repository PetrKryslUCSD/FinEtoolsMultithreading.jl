# This is about twenty percent faster than the original version.
function _populate_with_dofs!(dofs!, n, neighbors, dofnums, start)
    s1 = start[dofnums[n, 1]]
    p = 0
    @inbounds for k in neighbors
        for d in axes(dofnums, 2)
            dofs![s1+p] = dofnums[k, d]
            p += 1
        end
    end
    bl = p
    sort!(@view(dofs![s1:s1+bl-1]))
    @inbounds for d = 2:size(dofnums, 2)
        s = start[dofnums[n, d]]
        copy!(@view(dofs![s:s+bl-1]), @view(dofs![s1:s1+bl-1]))
    end
    return nothing
end

function _dof_block_lengths(IT, map, dofnums)
    nd = size(dofnums, 2)
    total_dofs = length(map) * nd
    lengths = _zeros_via_calloc(IT, total_dofs + 1)
    lengths[1] = 1
    @inbounds Threads.@threads for k in eachindex(map)
        kl = length(map[k]) * nd
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

function _prepare_arrays(IT, FT, map, dofnums)
    # First we create an array of the lengths of the dof blocks
    start = _dof_block_lengths(IT, map, dofnums)
    # Now we start overwriting the "lengths" array with the starts
    ThreadedScans.scan!(+, start) # equivalent to _acc_start_ptr!(start)
    sumlen = start[end] - 1
    dofs = Vector{IT}(undef, sumlen) # This will get filled in later
    nzval = _zeros_via_calloc(FT, sumlen) # This needs to be initialized for future accumulation
    return start, dofs, nzval
end

"""
    sparse_symmetric_csc_pattern(dofnums::Array{IT,2}, nrowscols, n2n, z = zero(Float64)) where {IT<:Integer}
    
Create symmetric sparse zero matrix (sparsity pattern).

Uses the following data structures:
```
    n2n = FENodeToNeighborsMap(n2e, fes.conn)
```
"""
function sparse_symmetric_csc_pattern(
    dofnums::Array{IT,2},
    nrowscols,
    n2n,
    z = zero(Float64),
) where {IT<:Integer}
    @assert length(n2n.map) == size(dofnums, 1)
    FT = typeof(z)
    # This is about an order of magnitude less expensive than the next step
    start, dofs, nzval = _prepare_arrays(IT, FT, n2n.map, dofnums)
    # This stops scaling for nthreads >= 32
    @inbounds Base.Threads.@threads for n in axes(dofnums, 1)
        _populate_with_dofs!(dofs, n, n2n.map[n], dofnums, start)
    end
    K = SparseMatrixCSC(nrowscols, nrowscols, start, dofs, nzval)
    return K
end
