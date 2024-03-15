using InteractiveUtils
# This is about twenty percent faster than the original version.
function _populate_arrays!(dofs!, nzval!, n, neighbors, dofnums, start)
    z = zero(eltype(nzval!))
    s1 = start[dofnums[n, 1]]
    p = 0
    @inbounds for k in neighbors
        for d in axes(dofnums, 2)
            dofs![s1+p] = dofnums[k, d]
            nzval![s1+p] = z
            p += 1
        end
    end
    bl = p
    sort!(@view(dofs![s1:s1+bl-1]))
    @inbounds for d in 2:size(dofnums, 2)
        s = start[dofnums[n, d]]
        for p in 0:1:bl-1
            dofs![s+p] = dofs![s1+p]
            nzval![s+p] = z
        end
    end
    return nothing
end

function _rowcol_lengths(IT, total_dofs, map, dofnums)
    nd = size(dofnums, 2)
    # lengths = Vector{IT}(undef, total_dofs + 1)
    lengths = fill(zero(IT), total_dofs + 1)
    @inbounds Threads.@threads for k in eachindex(map)
        kl = length(map[k]) * nd
        for d in axes(dofnums, 2)
            j = dofnums[k, d]
            lengths[j] = kl
        end
    end
    lengths[end] = 0
    return lengths
end

function _calculate_start(IT, map, dofnums)
    nd = size(dofnums, 2)
    total_dofs = length(map) * nd
    start = _rowcol_lengths(IT, total_dofs, map, dofnums)
    # Now we start overwriting the "lengths" array with the starts
    sumlen = 0
    len = start[1]
    sumlen += len
    start[1] = 1
    plen = len
    @inbounds for k in 2:total_dofs
        len = start[k]
        sumlen += len
        start[k] = start[k-1] + plen
        plen = len
    end
    start[end] = sumlen + 1
    return start
end

function _prepare_arrays(IT, FT, map, dofnums)
    # @code_warntype _calculate_start(IT, map, dofnums)
    @time start = _calculate_start(IT, map, dofnums)
    sumlen = start[end] - 1
    dofs = Vector{IT}(undef, sumlen)
    nzval = Vector{eltype(FT)}(undef, sumlen)
    return start, dofs, nzval
end

# function _prepare_arrays(IT, FT, n2n, dofnums)
#     nd = size(dofnums, 2)
#     total_dofs = length(n2n.map) * nd
#     lengths = Vector{IT}(undef, total_dofs + 1)
#     @time @inbounds for k in eachindex(n2n.map)
#         kl = length(n2n.map[k]) * nd
#         for d in axes(dofnums, 2)
#             j = dofnums[k, d]
#             lengths[j] = kl
#         end
#     end
#     lengths[end] = 0
#     # Now we start overwriting the lengths array with the starts
#     start = lengths
#     sumlen = 0
#     len = start[1]
#     sumlen += len
#     start[1] = 1
#     plen = len
#     @time @inbounds for k in 2:total_dofs
#         len = start[k]
#         sumlen += len
#         start[k] = start[k-1] + plen
#         plen = len
#     end
#     start[end] = sumlen + 1
#     dofs = Vector{IT}(undef, sumlen)
#     nzval = Vector{eltype(FT)}(undef, sumlen)
#     return start, dofs, nzval
# end

"""
    sparse_symmetric_zero(u, n2n, kind = :CSC)

Create symmetric sparse zero matrix (sparsity pattern).

Uses the following data structures:
```
    n2e = FENodeToFEMap(fes.conn, nnodes(u))
    n2n = FENodeToNeighborsMap(n2e, fes.conn)
```
"""
function sparse_symmetric_zero(u, n2n, kind = :CSC)
    @assert kind in [:CSC, :CSR]
    @assert length(n2n.map) == size(u.dofnums, 1)
    IT = eltype(u.dofnums)
    FT = eltype(u.values)
    nrowscols = nalldofs(u)
    start, dofs, nzval = _prepare_arrays(IT, FT, n2n.map, u.dofnums)
    @inbounds Base.Threads.@threads for n in axes(u.dofnums, 1)
        _populate_arrays!(dofs, nzval, n, n2n.map[n], u.dofnums, start)
    end
    if kind == :CSC
        K = _csc_matrix(start, dofs, nrowscols, nzval)
    elseif kind == :CSR
        K = _csr_matrix(start, dofs, nrowscols, nzval)
    end
    return K
end

function _csc_matrix(start, dofs, nrowscols, nzval)
    return SparseMatrixCSC(
        nrowscols,
        nrowscols,
        start,
        dofs,
        nzval,
    )
end

function _csr_matrix(start, dofs, nrowscols, nzval)
    return SparseMatricesCSR.SparseMatrixCSR{1}(
        nrowscols,
        nrowscols,
        start,
        dofs,
        nzval,
    )
end

