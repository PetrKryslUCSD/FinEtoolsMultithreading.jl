
function _populate_dofs(n, n2n, dofnums, start, dofs)
    nd = size(dofnums, 2)
    totd = length(n2n.map[n]) * nd
    _dofs = fill(zero(eltype(dofs)), totd)
    p = 1
    for k in n2n.map[n]
        for d in axes(dofnums, 2)
            _dofs[p] = dofnums[k, d]
            p += 1
        end
    end
    sort!(_dofs)
    for d in axes(dofnums, 2)
        j = dofnums[n, d]
        s = start[j]
        p = 0
        for m in eachindex(_dofs)
            dofs[s+p] = _dofs[m]
            p += 1
        end
    end
    return nothing
end

function _prepare_start_dofs(IT, n2n, dofnums)
    nd = size(dofnums, 2)
    total_dofs = length(n2n.map) * nd
    lengths = Vector{IT}(undef, total_dofs+1)
    for k in eachindex(n2n.map)
        kl = length(n2n.map[k]) * nd
        for d in axes(dofnums, 2)
            j = dofnums[k, d]
            lengths[j] = kl
        end
    end
    lengths[end] = 0
    # Now we start overwriting the lengths array with the starts
    start = lengths
    sumlen = 0
    len = start[1]
    sumlen += len
    start[1] = 1
    plen = len
    for k in 2:total_dofs
        len = start[k]
        sumlen += len
        start[k] = start[k-1] + plen 
        plen = len
    end
    start[end] = sumlen+1
    dofs = Vector{IT}(undef, sumlen)
    return start, dofs
end

"""
    sparsity_pattern_symmetric(fes, u)

Create symmetric sparsity pattern.
"""
function sparsity_pattern_symmetric(fes, u)
    IT = eltype(u.dofnums)
    n2e = FENodeToFEMap(fes.conn, nnodes(u))
    n2n = FENodeToNeighborsMap(n2e, fes.conn)
    start, dofs = _prepare_start_dofs(IT, n2n, u.dofnums)
    Base.Threads.@threads for n in axes(u.dofnums, 1)
        _populate_dofs(n, n2n, u.dofnums, start, dofs)
    end
    return start, dofs
end

function csc_matrix_pattern(fes, u)
    start, dofs = sparsity_pattern_symmetric(fes, u)
    return SparseMatrixCSC(nalldofs(u), nalldofs(u), start, dofs, fill(zero(eltype(u.values)), length(dofs)))
end