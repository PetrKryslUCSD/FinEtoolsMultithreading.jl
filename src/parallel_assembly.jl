using FinEtools
using SparseArrays
import FinEtools.AssemblyModule.eltype
import FinEtools.AssemblyModule.startassembly!
import FinEtools.AssemblyModule.assemble!
import FinEtools.AssemblyModule.makematrix!

"""
    SysmatAssemblerSparsePatt{IT, MBT, IBT} <: AbstractSysmatAssembler

Type for assembling a sparse global matrix from elementwise matrices using a
sparsity pattern (sparse matrix with all potential non-zeroes represented, but
all those numbers are actually zero).

!!! note

    All fields of the datatype are private. The type is manipulated by the
    functions `startassembly!`, `assemble!`, and `makematrix!`.
"""
mutable struct SysmatAssemblerSparsePatt{TPATT<:SparseMatrixCSC} <: AbstractSysmatAssembler
    _pattern::TPATT
    _nomatrixresult::Bool
    _force_init::Bool
end

function associate_pattern(self::T, pattern::TPATT) where {T<:SysmatAssemblerSparsePatt,TPATT<:SparseMatrixCSC}
    self._pattern = pattern
end

"""
    SysmatAssemblerSparsePatt(z = zero(T)) where {T}

Construct a sparse system matrix assembler.

The matrix entries are of type `T`. 

# Example

This is how a sparse matrix is assembled from two rectangular dense matrices.
```
    a = SysmatAssemblerSparsePatt(0.0)
    associate_pattern(a, sparse(zeros(7, 7)))
    startassembly!(a, 5, 5, 3, 7, 7)
    m = [0.24406   0.599773    0.833404  0.0420141
        0.786024  0.00206713  0.995379  0.780298
        0.845816  0.198459    0.355149  0.224996]
    assemble!(a, m, [1 7 5], [5 2 1 4])
    m = [0.146618  0.53471   0.614342    0.737833
         0.479719  0.41354   0.00760941  0.836455
         0.254868  0.476189  0.460794    0.00919633
         0.159064  0.261821  0.317078    0.77646
         0.643538  0.429817  0.59788     0.958909]
    assemble!(a, m, [2 3 1 7 5], [6 7 3 4])
    A = makematrix!(a)
```
Here `A` is a sparse matrix of the size 7x7.

!!! note

    The `nomatrixresult` field is ignored.

"""
function SysmatAssemblerSparsePatt(z::T) where {T}
    return SysmatAssemblerSparsePatt(
        spzeros(T, 1, 1),
        false,
        false,
    )
end

function SysmatAssemblerSparsePatt()
    return SysmatAssemblerSparsePatt(zero(Float64))
end

function eltype(self::TPATT) where {TPATT<:SysmatAssemblerSparsePatt}
    eltype(self._pattern.nzval)
end

"""
    startassembly!(self::SysmatAssemblerSparsePatt{T},
        elem_mat_nrows::IT,
        elem_mat_ncols::IT,
        n_elem_mats::IT,
        row_nalldofs::IT,
        col_nalldofs::IT;
        force_init = false
        ) where {T, IT<:Integer}

Start the assembly of a global matrix.

The method makes buffers for matrix assembly. It must be called before
the first call to the method `assemble!`.

# Arguments
- `elem_mat_nrows` = row dimension of the element matrix;
- `elem_mat_ncols` = column dimension of the element matrix;
- `n_elem_mats` = number of element matrices;
- `row_nalldofs`= The total number of rows as a tuple;
- `col_nalldofs`= The total number of columns as a tuple.

This is a noop. The pattern has already been built, per our assumption.

# Returns
- `self`: the modified assembler.

"""
function startassembly!(
    self::TPATT,
    elem_mat_nrows::IT,
    elem_mat_ncols::IT,
    n_elem_mats::IT,
    row_nalldofs::IT,
    col_nalldofs::IT;
    force_init = false,
) where {TPATT<:SysmatAssemblerSparsePatt,IT<:Integer}
    return self
end

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

function _updroworcol!(nzval, i, v, st, fi, r_or_c)
    k = _binary_search(r_or_c, i, st, fi)
    if k > 0
        nzval[k] += v
    end
end

"""
    assemble!(
        self::SysmatAssemblerSparsePatt,
        mat::MT,
        dofnums_row::IT,
        dofnums_col::IT,
    ) where {MT, IT}

Assemble a rectangular matrix.
"""
function assemble!(
    self::TPATT,
    mat::MBT,
    dofnums_row::CIT,
    dofnums_col::CIT,
) where {TPATT<:SysmatAssemblerSparsePatt,MBT,CIT}
    # Assembly of a rectangular matrix.
    # The method assembles a rectangular matrix using the two vectors of
    # equation numbers for the rows and columns.
    nrows = length(dofnums_row)
    ncolumns = length(dofnums_col)
    size(mat) == (nrows, ncolumns) || error("Wrong size of matrix")
    row_nalldofs, col_nalldofs = size(self._pattern)
    for j = 1:ncolumns
        dj = dofnums_col[j]
        dj < 1 && error("Column degree of freedom < 1")
        dj > col_nalldofs && error("Column degree of freedom > size")
        for i = 1:nrows
            di = dofnums_row[i]
            di < 1 && error("Row degree of freedom < 1")
            di > row_nalldofs && error("Row degree of freedom > size")
            _updroworcol!(
                self._pattern.nzval,
                di,
                mat[i, j],
                self._pattern.colptr[j],
                self._pattern.colptr[j+1] - 1,
                self._pattern.rowval,
            )
        end
    end
    return self
end

"""
    makematrix!(self::SysmatAssemblerSparsePatt)

Make a sparse matrix.

The sparse matrix (i.e. a sparsity pattern with the nonzero values filled in) is returned.
 
"""
function makematrix!(self::TPATT) where {TPATT<:SysmatAssemblerSparsePatt}
    return self._pattern
end

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
    start_assembler!(
        assembler::AT
    ) where {AT<:AbstractSysmatAssembler}

Start parallel assembly on the global assembler.

The assembler is informed that no matrix result should be computed. All the data
should be preserved in the assembly buffers.
"""
function start_parallel_assembler!(assembler::AT, force_init = false) where {AT<:SysmatAssemblerSparsePatt}
    if force_init
        fill!(zero(eltype(assembler._pattern.nzval)), assembler._pattern.nzval)
    end
    return assembler
end

"""
    parallel_matrix_assembly(
        femms::AbstractArray{FEMM, 1},
        assemblers::AbstractVector{AT},
        matrixcomputation!::F,
        kind = :threaded,
    ) where {FEMM<:AbstractFEMM, AT<:AbstractSysmatAssembler, F<:Function}

Execute the assembly in parallel.
"""
function parallel_matrix_assembly!(
    assembler,
    decomposition,
    matrixcomputation!::F,
    ntasks = Threads.nthreads()
) where {F<:Function}
    for d in decomposition
        femms = d
        Threads.@threads for j in eachindex(femms)
            matrixcomputation!(femms[j], assembler)
        end
    end
    return true
end
