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

function eltype(self::TPATT) where {TPATT<:SysmatAssemblerSparsePatt}
    eltype(self._pattern.nzval)
end

function SysmatAssemblerSparsePatt(pattern::TPATT) where {TPATT<:SparseMatrixCSC}
    return SysmatAssemblerSparsePatt(pattern, false, false)
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
    if force_init
        fill!(zero(eltype(self._pattern.nzval)), self._pattern.nzval)
    end
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
                self._pattern.colptr[dj],
                self._pattern.colptr[dj+1] - 1,
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
