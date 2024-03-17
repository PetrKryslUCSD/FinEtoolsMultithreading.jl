
"""
    subdomainfemms(fes, ntasks, crsubdom)

Create finite element machines for the subdomains.

# Arguments
- `fes` = finite element set 
- `ntasks` = number of tasks (subdomains)
- `crsubdom` = function to create one finite element machine
    Example: 
    ```
    function createsubdomain(fessubset)
        FEMMDeforLinear(MR, IntegDomain(fessubset, GaussRule(3, 2)), material)
    end
    ```
"""
function subdomainfemms(fes, ntasks, crsubdom)
    if ntasks == 1
        return [crsubdom(fes)]
    else
        sds = []
        chks = chunks(1:count(fes), ntasks)
        resize!(sds, ntasks)
        Threads.@threads for k = 1:ntasks
            (ch, j) = chks[k]
            sds[k] = crsubdom(subset(fes, ch))
        end
        return typeof(sds[1])[sds[k] for k in eachindex(sds)]
    end
end
