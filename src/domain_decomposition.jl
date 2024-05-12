"""
    decompose(fes, coloring, createsubd, ntasks = Threads.nthreads())

Create finite element machines for the subdomains.

# Arguments
- `fes` = finite element set,
- `ntasks` = number of tasks (subdomains),
- `element_colors` = array of element colors, one for each element,
- `unique_colors` = array of the unique element colours,
- `createsubd` = function to create one finite element machine 
    
    Example: 
    ```
    function createsubd(fessubset)
        FEMMDeforLinear(MR, IntegDomain(fessubset, GaussRule(3, 2)), material)
    end
    ```

Create a vector of vectors: for each unique color in the element coloring, the
vector of the finite element machines is stored. The finite element machines are
created using the function `createsubd`.

The matrix is created by going sequentially through the unique colors
and then in parallel execute all the finite element machines for that color.
"""
function decompose(fes, coloring, createsubd,
    ntasks=Threads.nthreads())
    el_colors, uniq_colors = coloring
    decomp = fill([], length(uniq_colors))
    Threads.@threads for i in eachindex(uniq_colors)
        c = uniq_colors[i]
        ellist = findall(_c -> _c == c, el_colors)
        _fes = subset(fes, ellist)
        decomp[i] = _make_femms(_fes, ntasks, createsubd)
    end
    return decomp
end

function _make_femms(fesofacolor, ntasks, createsubd)
    chks = chunks(1:count(fesofacolor), ntasks)
    return [createsubd(subset(fesofacolor, ch)) for (ch, j) in chks]
end

# using LazyArrays

# i = 1000
# is = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]

# function f(i, is)
#     for n in ApplyArray(vcat, i, is)
#         n
#     end
# end

# @time f(i, is)