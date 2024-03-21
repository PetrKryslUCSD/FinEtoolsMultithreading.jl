
"""
    subdomainfemms(fes, ntasks, element_colors, unique_colors, crsubdom)

Create finite element machines for the subdomains.

# Arguments
- `fes` = finite element set,
- `ntasks` = number of tasks (subdomains),
- `element_colors` = array of element colors, one for each element,
- `unique_colors` = array of the unique element colours,
- `crsubdom` = function to create one finite element machine
    Example: 
    ```
    function createsubdomain(fessubset)
        FEMMDeforLinear(MR, IntegDomain(fessubset, GaussRule(3, 2)), material)
    end
    ```
"""
function domain_decomposition(fes, ntasks, element_colors, unique_colors, crsubdom)
    decomposition = []
    Threads.@threads for color in element_colors
        color_list = findall(x -> x == color, element_colors)
        subsetfes = subset(fes, color_list)
        chks = chunks(1:count(subsetfes), ntasks)
        sds = []
        resize!(sds, ntasks)
        for k in eachindex(chks)
            (ch, j) = chks[k]
            sds[k] = crsubdom(subset(fes, ch))
        end
        push!(decomposition, typeof(sds[1])[sds[k] for k in eachindex(sds)])
    end
    return decomposition
end
