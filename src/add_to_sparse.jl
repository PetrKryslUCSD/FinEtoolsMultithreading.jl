"""
    zerooutsparse(S)

Zero out the stored entries of the matrix. The sparsity pattern is not affected. 
"""
function zerooutsparse(S)
    S.nzval .= zero(eltype(S.nzval))
    return S
end

function _binary_search(array::Array{IT,1}, target::IT, left::IT, right::IT)  where {IT}
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

function _updcol!(nzval, i, v, st, fi, rowval)
    k = _binary_search(rowval, i, st, fi)
    if k > 0
        nzval[k] += v
    end
end

"""
    addtosparse(S, I, J, V)

Add the values from the array `V` given the row and column indexes in the arrays
`I` and `J`. The expectation is that the indexes respect the sparsity pattern of
the sparse array `S`. 
"""
function addtosparse(S, I, J, V)
    nzval = S.nzval; colptr = S.colptr; rowval = S.rowval;
    Threads.@threads for t in eachindex(J)
        j = J[t]
        _updcol!(nzval, I[t], V[t], colptr[j], colptr[j+1]-1, rowval)
    end
    return S
end
