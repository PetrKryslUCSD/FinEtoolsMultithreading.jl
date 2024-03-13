
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
    return [crsubdom(subset(fes, ch)) for (ch, j) in chunks(1:count(fes), ntasks)]
end
