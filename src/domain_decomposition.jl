"""
    domain_decomposition(fes, coloring, createsubd, ntasks = Threads.nthreads())

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
function domain_decomposition(fes, coloring, createsubd,
    ntasks=Threads.nthreads())
    el_colors, uniq_colors = coloring
    decomposition = fill([], length(uniq_colors))
    Threads.@threads for i in eachindex(uniq_colors)
        ellist = findall(c -> c == uniq_colors[i], el_colors)
        fesofacolor = subset(fes, ellist)
        decomposition[i] =
            _make_femms(fesofacolor, ntasks, createsubd)
    end
    return decomposition
end

function _make_femms(fesofacolor, ntasks, createsubd)
    chks = chunks(1:count(fesofacolor), ntasks)
    return [createsubd(subset(fesofacolor, ch)) for (ch, j) in chks]
end
