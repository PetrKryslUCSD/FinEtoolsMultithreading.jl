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
function domain_decomposition(fes, coloring, createsubd, ntasks=Threads.nthreads())
    element_colors, unique_colors = coloring
    decomposition = []
    resize!(decomposition, length(unique_colors))
    Threads.@threads for i in eachindex(unique_colors)
        decomposition[i] =
            _make_femms(fes, ntasks, createsubd,
                findall(c -> c == unique_colors[i], element_colors)
            )
    end
    return decomposition
end

function _make_femms(fes, ntasks, createsubd, color_list)
    subsetfes = subset(fes, color_list)
    chks = chunks(1:count(subsetfes), ntasks)
    sds = []
    resize!(sds, length(chks))
    for k = 1:length(chks)
        (ch, j) = chks[k]
        sds[k] = createsubd(subset(subsetfes, ch))
    end
    return typeof(sds[1])[sds[k] for k in eachindex(sds)]
end
