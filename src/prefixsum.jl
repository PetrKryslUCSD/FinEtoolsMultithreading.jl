using ThreadedScans
using ChunkSplitters


function _scan!(arr)
    n = length(arr)
    @inbounds for i in 2:n
        arr[i] += arr[i-1]
    end
    return arr
end

"""
    psp_scan!(arr)

Perform prefix sum scan on the input array `arr` in place.

# Arguments
- `arr`: The input array to perform prefix sum scan on. The array is modified
  in-place.

The computation is multithreaded, and uses all available threads.
"""
function psp_scan!(arr)
    n = length(arr)
    ntasks = Threads.nthreads()
    chks = chunks(1:n, ntasks)
    Threads.@threads for s in 1:length(chks)
        _scan!(@view arr[chks[s][1]])
    end
    segment_ends = [arr[chks[s][1][end]] for s in 1:length(chks)]
    _scan!(segment_ends)
    Threads.@threads for s in 2:length(chks)
        @view(arr[chks[s][1]]) .+= segment_ends[s-1]
    end
    return arr
end
