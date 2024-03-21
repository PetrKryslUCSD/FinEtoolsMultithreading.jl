module PQuickSort

using LinearAlgebra

@inline function _partition!(A, perm, pivot, left,right)
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

function __pquicksortperm!(A, perm, i, j, task)
    if task <= 0
        quicksortperm!(A, perm, i, j)
        return
    end
    left, right = _partition!(A, perm, A[(j+i) >>> 1], i, j)
    t = Threads.@spawn __pquicksortperm!(A, perm, $i, $right, task - 1)
    __pquicksortperm!(A, perm, left, j, task - 1)
    wait(t)
    return
end

# This has the same order of the arguments as the built in sortperm! 
function pquicksortperm!(perm, A, ntasks = Threads.nthreads())
    __pquicksortperm!(A, perm, 1, length(A), ntasks)
end

end # module

# using LinearAlgebra
# using Main.PQuickSort
# function g()
#     N = 100000000
#     # N = 10
#     a = collect((N):-1:1)
#     p = similar(a)
#     @time sortperm!(p, a)
#     # @show a, prm
#     aref = sort(deepcopy(a))

#     a = collect((N):-1:1)
#     prm = collect(1:length(a))
#     @time PQuickSort.pquicksortperm!(prm, a)
#     # @show a, prm

#     norm(a - sort(aref)) == 0

    

#     a = collect((N):-1:1)
#     prm = collect(1:length(a))
#     @time PQuickSort.pquicksortperm!(prm, a)
# end