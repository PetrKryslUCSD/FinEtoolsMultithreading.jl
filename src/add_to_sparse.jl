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

"""
    addtosparse(S::T, I, J, V) where {T<:SparseArrays.SparseMatrixCSC}

Add values to sparse CSC matrix.

Add the values from the array `V` given the row and column indexes in the arrays
`I` and `J`. The expectation is that the indexes respect the sparsity pattern of
the sparse array `S`. 
"""
function addtosparse(S::T, I, J, V, ntasks) where {T<:SparseArrays.SparseMatrixCSC}
    nzval = S.nzval
    colptr = S.colptr
    rowval = S.rowval
    blockl = min(length(I), 7)
    Threads.@sync begin
        for t in 1:ntasks
            Threads.@spawn let task = $t
                @show task
                for b in 1:blockl
                    for s in b:blockl:length(J)
                        j = J[s]
                        if rem(j, ntasks) + 1 == task
                            _updroworcol!(nzval, I[s], V[s], colptr[j], colptr[j+1] - 1, rowval)
                        end
                    end
                end
            end
        end
    end
    
    return S
end

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
    add_to_matrix!(
        S,
        assembler::AT
    ) where {AT<:AbstractSysmatAssembler}

Update the global matrix.

Use the sparsity pattern in `S`, and the COO data collected in the assembler.
"""
function add_to_matrix!(S, assembler::AT; ntasks = Threads.nthreads()) where {AT<:AbstractSysmatAssembler}
    # At this point all data is in the buffer
    assembler._buffer_pointer = assembler._buffer_length + 1
    setnomatrixresult(assembler, false)
    return addtosparse(S, assembler._rowbuffer, assembler._colbuffer, assembler._matbuffer, ntasks)
end

function do_one_entry(c)
    while true
        k = take!(c)
        if k == 0
            break
        end
        j = J[k]
        _updroworcol!(nzval, I[k], V[k], colptr[j], colptr[j+1] - 1, rowval)
    end
end

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
