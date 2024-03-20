module PQuickSort

const SEQ_THRESH = 2^17

@noinline function _partition!(A, perm, pivot, left,right)
    @inbounds while left <= right
        while A[left] < pivot
            left += 1
        end
        while A[right] > pivot
            right -= 1
        end
        if left <= right
            A[left], A[right] = A[right], A[left]
            perm[left], perm[right] = perm[right], perm[left]
            left += 1
            right -= 1
        end
    end
    return (left,right)
end

function quicksortperm!(A, perm, i=1, j=length(A))
    if j > i
        left, right = _partition!(A, perm, A[(j+i) >>> 1], i, j)
        quicksortperm!(A, perm, i, right)
        quicksortperm!(A, perm, left, j)
    end
end

function pquicksortperm!(A, perm, i=1, j=length(A))
    if j-i <= SEQ_THRESH
        quicksortperm!(A, perm, i, j)
        return
    end
    left, right = _partition!(A, perm, A[(j+i) >>> 1], i, j)
    # t = Threads.@spawn pquicksortperm!(A, perm, $i, $right)
    # pquicksortperm!(A, perm, left, j)
    # wait(t)
    Threads.@threads for w in  ((i, right), (left, j))
        pquicksortperm!(A, perm, w[1], w[2])
    end
    return
end

end # module
