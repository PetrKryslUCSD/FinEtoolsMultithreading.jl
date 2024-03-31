module m

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

N = 10000
r = collect(1:N)
g = collect(N:-1:1)
NLOOP = 1000
d = Dict{Int, Int}(zip(r, r))

function test_bsearch(N, r, g)
    for l in 1:NLOOP
        for R in Int.(round.([0.13, 0.39, 0.77, 0.98] * N))
            v = do_bsearch(N, r, g, R)
            @assert v == g[R]
        end
    end
    nothing
end

function do_bsearch(N, r, g, R)
    k = _binary_search(r, R, 1, N)
    return g[k]
end

function test_dict(N, d, g)
    for l in 1:NLOOP
        for R in Int.(round.([0.13, 0.39, 0.77, 0.98] * N))
            v = do_dict(N, d, g, R)
            @assert v == g[R]
        end
    end
    nothing
end

function do_dict(N, d, g, R)
    k = d[R]
    return g[k]
end

using BenchmarkTools

@btime test_bsearch(N, r, g)
@btime test_bsearch($N, $r, $g)

@btime test_dict(N, d, g)
@btime test_dict($N, $d, $g)

end  # module
