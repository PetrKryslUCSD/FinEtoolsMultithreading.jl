using FinEtools
using SparseArrays
import FinEtools.AssemblyModule.eltype
import FinEtools.AssemblyModule.startassembly!
import FinEtools.AssemblyModule.assemble!
import FinEtools.AssemblyModule.makematrix!

include("assembler_patt.jl")
include("assembler_patt_w_lookup.jl")

function _check_femm_compatibility(femms::AbstractArray{FEMM,1}) where {FEMM<:AbstractFEMM}
    for j in eachindex(femms)
        iselementbased(femms[j]) || error("FEMM is not element-based")
        (nameof(typeof(femms[j])) === nameof(typeof(femms[1]))) ||
            error("All FEMMs must be of the same type")
        (
            nameof(typeof(finite_elements(femms[j]))) ===
            nameof(typeof(finite_elements(femms[1])))
        ) || error("All finite elements must be of the same type")
    end
    return true
end

"""
    parallel_matrix_assembly!(
        assmblr,
        decomposition,
        matrixupdt!::F
    ) where {F<:Function}

Execute the assembly in parallel.

The decomposition is a vector of vector of FEMMs.
As many tasks as there are FEMMs at any level are spawned.

The function `matrixupdt!` updates the assembler.
"""
function parallel_matrix_assembly!(
    assmblr,
    decomposition,
    matrixupdt!::F,
) where {F<:Function}
    for femms in decomposition
        Threads.@sync begin
            for femm in femms
                Threads.@spawn matrixupdt!(femm, assmblr)
            end
        end
    end
    return assmblr._pattern
end
