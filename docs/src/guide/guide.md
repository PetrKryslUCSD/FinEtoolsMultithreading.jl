# Guide

The [`FinEtools`](https://petrkryslucsd.github.io/FinEtools.jl/latest/index.html) 
package is used here to solve a variety of finite element problems.
This package can provide an overlay to parallelize sparse matrix assembly.

There is one high-level function that can be used to parallelize the assembly of any sparse matrix -- acoustic mass or stiffness,  conductivity matrix, stiffness or mass matrix, etc.

This bit of code will assemble the diffusion by linear forms stiffness matrix stored in the CSR format:
```
function createsubdomain(fessubset)
    FEMMBase(IntegDomain(fessubset, GaussRule(3, 2)))
end

function matrixcomputation!(femm, assembler)
    bilform_diffusion(femm, assembler, geom, psi, DataCache(Matrix(1.0 * LinearAlgebra.I(3))))
end

K1 = parallel_make_matrix(fes, psi, createsubdomain, matrixcomputation!;
    ntasks=ntasks, kind=:CSR)
```
The user can ask for any number of tasks to be used (even though it would be best to match it to the number of available threads).