"""
    domain_decomposition(fes, ntasks, element_colors, unique_colors, createsubdomain)

Create finite element machines for the subdomains.

# Arguments
- `fes` = finite element set,
- `ntasks` = number of tasks (subdomains),
- `element_colors` = array of element colors, one for each element,
- `unique_colors` = array of the unique element colours,
- `createsubdomain` = function to create one finite element machine 
    
    Example: 
    ```
    function createsubdomain(fessubset)
        FEMMDeforLinear(MR, IntegDomain(fessubset, GaussRule(3, 2)), material)
    end
    ```

Create a vector of vectors: for each unique color in the element coloring, the
vector of the finite element machines is stored. The finite element machines are
created using the function `createsubdomain`.

The matrix is created by going sequentially through the unique colors
and then in parallel execute all the finite element machines for that color.
"""
function domain_decomposition(fes, ntasks, element_colors, unique_colors, createsubdomain)
    decomposition = []
    Threads.@threads for color in unique_colors
        color_list = findall(x -> x == color, element_colors)
        subsetfes = subset(fes, color_list)
        chks = chunks(1:count(subsetfes), ntasks)
        sds = []
        resize!(sds, length(chks))
        for k in 1:length(chks)
            (ch, j) = chks[k]
            sds[k] = createsubdomain(subset(subsetfes, ch))
        end
        push!(decomposition, typeof(sds[1])[sds[k] for k in eachindex(sds)])
    end
    return decomposition
end
