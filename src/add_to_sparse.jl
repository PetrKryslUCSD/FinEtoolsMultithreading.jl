"""
    zerooutsparse(S)

Zero out the stored entries of the matrix. The sparsity pattern is not affected. 
"""
function zerooutsparse(S)
    S.nzval .= zero(eltype(S.nzval))
    return S
end

function _binary_search(array::Array{IT,1}, target::IT, left::IT, right::IT) where {IT}
    @inbounds while left <= right # Generating the middle element position 
        mid = fld((left + right), 2) # If element > mid, then it can only be present in right subarray
        if array[mid] < target
            left = mid + 1 # If element < mid, then it can only be present in left subarray 
        elseif array[mid] > target
            right = mid - 1 # If element is present at the middle itself 
        else # == 
            return mid
        end
    end
    return 0
end

function _updroworcol!(nzval, i, v, st, fi, r_or_c)
    k = _binary_search(r_or_c, i, st, fi)
    if k > 0
        nzval[k] += v
    end
end

function _find_breaks(J, ntasks)
    chks = chunks(1:length(J), ntasks)
    lo = fill(1, ntasks)
    hi = fill(length(J), ntasks)
    for t in 1:ntasks
        ch = chks[t]
        rowfirst = minimum(ch[1])
        rowlast = maximum(ch[1])
        lo[t] = rowfirst
        hi[t] = rowlast
    end
    for t in 2:ntasks
        p = hi[t-1]
        c = J[p]
        while (J[p] == c) p += 1; end 
        hi[t-1] = p - 1
        lo[t] = p
    end
    return lo, hi
end

# function addtosparse(S::T, I, J, V, ntasks) where {T<:SparseArrays.SparseMatrixCSC}
#     nzval = S.nzval
#     colptr = S.colptr
#     rowval = S.rowval
#     @time prm = sortperm(J)
#     @time begin
#     I = I[prm]
#     J = J[prm]
#     V = V[prm]
#     end
#     lo, hi = _find_breaks(J, ntasks)
#     # for k in 1:length(lo)
#     #     @show J[prm[lo[k]]], J[prm[max(1,lo[k]-3):lo[k]+3]]    
#     #     @show J[prm[hi[k]]], J[prm[hi[k]-3:min(length(J),hi[k]+3)]]    
#     # end
#     # Threads.@sync begin
#     #     for t in 1:ntasks
#     #        Threads.@spawn let l = lo[$t], h = hi[$t]
#     #         @show l:h
#     #             @inbounds for s in l:h
#     #                 j = J[s]
#     #                 _updroworcol!(nzval, I[s], V[s], colptr[j], colptr[j+1] - 1, rowval)
#     #             end
#     #         end
#     #     end
#     # end
#     Threads.@threads for t in 1:ntasks
#         @inbounds for s in lo[t]:hi[t]
#             j = J[s]
#             _updroworcol!(nzval, I[s], V[s], colptr[j], colptr[j+1] - 1, rowval)
#         end
#     end
#     return S
# end

function addtosparse(S::T, I, J, V, ntasks) where {T<:SparseArrays.SparseMatrixCSC}
    nzval = S.nzval
    colptr = S.colptr
    rowval = S.rowval
    @time prm = sortperm(J)
    @time begin
    I = I[prm]
    J = J[prm]
    V = V[prm]
    end
    lo, hi = _find_breaks(J, ntasks)
    # Threads.@threads for t in 1:ntasks
    #     @inbounds for s in lo[t]:hi[t]
    #         j = J[s]
    #         _updroworcol!(nzval, I[s], V[s], colptr[j], colptr[j+1] - 1, rowval)
    #     end
    # end
    Threads.@sync begin
        for t in 1:ntasks
            Threads.@spawn let 
                for s in lo[t]:hi[t]
                    j = J[s]
                    _updroworcol!(nzval, I[s], V[s], colptr[j], colptr[j+1] - 1, rowval)
                end
            end
        end
    end
    return S
end

"""
    addtosparse(S::T, I, J, V) where {T<:SparseArrays.SparseMatrixCSC}

Add values to sparse CSC matrix.

Add the values from the array `V` given the row and column indexes in the arrays
`I` and `J`. The expectation is that the indexes respect the sparsity pattern of
the sparse array `S`. 
"""
# function addtosparse(S::T, I, J, V, ntasks) where {T<:SparseArrays.SparseMatrixCSC}
#     nzval = S.nzval
#     colptr = S.colptr
#     rowval = S.rowval
#     chks = chunks(1:size(S, 2), ntasks)
#     Threads.@sync begin
#         for t in 1:ntasks
#             ch = chks[t]
#             rowfirst = minimum(ch[1])
#             rowlast = maximum(ch[1])
#             Threads.@spawn let rowfirst = $rowfirst, rowlast = $rowlast
#                 for s in eachindex(J)
#                     j = J[s]
#                     if (rowfirst <= j <= rowlast)
#                         _updroworcol!(nzval, I[s], V[s], colptr[j], colptr[j+1] - 1, rowval)
#                     end
#                 end
#             end
#         end
#     end
#     return S
# end

"""
    addtosparse(S::T, I, J, V) where {T<:SparseMatricesCSR.SparseMatrixCSR}

Add values to sparse CSR matrix.

Add the values from the array `V` given the row and column indexes in the arrays
`I` and `J`. The expectation is that the indexes respect the sparsity pattern of
the sparse array `S`. 
"""
function addtosparse(S::T, I, J, V, ntasks) where {T<:SparseMatricesCSR.SparseMatrixCSR}
    nzval = S.nzval
    rowptr = S.rowptr
    colval = S.colval
    chks = chunks(1:size(S, 1), ntasks)
    Threads.@sync begin
        for ch in chks
            from = minimum(ch[1])
            to = maximum(ch[1])
            Threads.@spawn let from = $from, to = $to
                @inbounds for t in eachindex(I)
                    i = I[t]
                    if (from <= i <= to) 
                        _updroworcol!(nzval, J[t], V[t], rowptr[i], rowptr[i+1] - 1, colval)
                    end
                end
            end
        end
    end
    return S
end

"""
    add_to_matrix!(S, assembler::AT, ntasks=Threads.nthreads()) where {AT<:AbstractSysmatAssembler}

Update the global matrix.

Use the sparsity pattern in `S`, and the COO data collected in the assembler.
"""
function add_to_matrix!(S, assembler::AT, ntasks=Threads.nthreads()) where {AT<:AbstractSysmatAssembler}
    # At this point all data is in the buffer
    assembler._buffer_pointer = assembler._buffer_length + 1
    setnomatrixresult(assembler, false)
    return addtosparse(S, assembler._rowbuffer, assembler._colbuffer, assembler._matbuffer, ntasks)
end

# function do_one_entry(c)
#     while true
#         k = take!(c)
#         if k == 0
#             break
#         end
#         j = J[k]
#         _updroworcol!(nzval, I[k], V[k], colptr[j], colptr[j+1] - 1, rowval)
#     end
# end

# ntasks = Threads.nthreads()
# c = [Channel{Int}() for _  in 1:ntasks]
# tasks = Task[]
# for t in 1:ntasks
#     push!(tasks, @task(do_one_entry(c[t])))
#     schedule(tasks[t])
# end
# for k in eachindex(J)
#     j = J[k]
#     t = rem(j, ntasks) + 1
#     put!(c[t], k)
# end
# for t in 1:ntasks
#     put!(c[t], 0)
# end

# sweepfrom = fill(typemax(eltype(I)), ntasks)
#     sweepto = fill(typemin(eltype(I)), ntasks)
#     @show J[1:20]
#     Threads.@threads for s in eachindex(J)
#         @inbounds for t in 1:ntasks
#             ch = chks[t]
#             rowfirst = minimum(ch[1])
#             rowlast = maximum(ch[1])
#             if J[s] >= rowfirst && J[s] <= rowlast
#                 sweepfrom[t] = min(sweepfrom[t], s)
#                 sweepto[t] = max(sweepto[t], s)
#             end
#         end
#     end
#     @show length(J), sweepfrom, sweepto


# function addtosparse(S::T, I, J, V, ntasks) where {T<:SparseArrays.SparseMatrixCSC}
#     nzval = S.nzval
#     colptr = S.colptr
#     rowval = S.rowval
#     IT = eltype(I)
#     blocksize = 100000
#     d = Dict{IT, Vector{Tuple{IT, IT}}}()
#     p = 0
#     while p < length(J)
#         ub = min(blocksize, length(J) - p)
#         for b in 1:ub
#             p += 1
#             j = J[p]
#             if !(j in keys(d))
#                 d[j] = Tuple{IT, IT}[]
#             end
#             push!(d[j], (p, I[p]))
#         end
#         Threads.@threads for j in [j for j in keys(d)]
#             for (_p, _i) in d[j]
#                 _updroworcol!(nzval, _i, V[_p], colptr[j], colptr[j+1] - 1, rowval)
#             end
#         end
#         empty!(d)
#     end
#     return S
# end
