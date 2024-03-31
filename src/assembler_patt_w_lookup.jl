using FinEtools
using SparseArrays
import FinEtools.AssemblyModule.eltype
import FinEtools.AssemblyModule.startassembly!
import FinEtools.AssemblyModule.assemble!
import FinEtools.AssemblyModule.makematrix!

"""
    SysmatAssemblerSparsePattwLookup{IT, MBT, IBT} <: AbstractSysmatAssembler

Type for assembling a sparse global matrix from elementwise matrices using a
sparsity pattern (sparse matrix with all potential non-zeroes represented, but
all those numbers are actually zero).

!!! note

    All fields of the datatype are private. The type is manipulated by the
    functions `startassembly!`, `assemble!`, and `makematrix!`.
"""
mutable struct SysmatAssemblerSparsePattwLookup{TPATT<:SparseMatrixCSC} <: AbstractSysmatAssembler
    _pattern::TPATT
    _rows::Vector{Dict{Int, Int}}
    _nomatrixresult::Bool
    _force_init::Bool
end

function eltype(self::TPATT) where {TPATT<:SysmatAssemblerSparsePattwLookup}
    eltype(self._pattern.nzval)
end

function SysmatAssemblerSparsePattwLookup(pattern::TPATT) where {TPATT<:SparseMatrixCSC}
    rows = Dict{Int, Int}[]
    resize!(rows, size(pattern, 2))
    Threads.@threads for c in axes(pattern, 2)
        rows[c] = Dict{Int, Int}()
        for k in pattern.colptr[c]:(pattern.colptr[c+1] - 1)
            r = pattern.rowval[k]
            rows[c][r] = k
        end
    end
    return SysmatAssemblerSparsePattwLookup(pattern, rows, false, false)
end

"""
    startassembly!(self::SysmatAssemblerSparsePattwLookup{T},
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
) where {TPATT<:SysmatAssemblerSparsePattwLookup,IT<:Integer}
    if force_init
        fill!(zero(eltype(self._pattern.nzval)), self._pattern.nzval)
    end
    return self
end

"""
    assemble!(
        self::SysmatAssemblerSparsePattwLookup,
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
) where {TPATT<:SysmatAssemblerSparsePattwLookup,MBT,CIT}
    # Assembly of a rectangular matrix.
    # The method assembles a rectangular matrix using the two vectors of
    # equation numbers for the rows and columns.
    nrows = length(dofnums_row)
    ncolumns = length(dofnums_col)
    size(mat) == (nrows, ncolumns) || error("Wrong size of matrix")
    row_nalldofs, col_nalldofs = size(self._pattern)
    @inbounds for j = 1:ncolumns
        dj = dofnums_col[j]
        # dj < 1 && error("Column degree of freedom < 1")
        # dj > col_nalldofs && error("Column degree of freedom > size")
        for i = 1:nrows
            di = dofnums_row[i]
            # di < 1 && error("Row degree of freedom < 1")
            # di > row_nalldofs && error("Row degree of freedom > size")
            k = self._rows[dj][di]
            self._pattern.nzval[k] += mat[i, j]
        end
    end
    return self
end

"""
    makematrix!(self::SysmatAssemblerSparsePattwLookup)

Make a sparse matrix.

The sparse matrix (i.e. a sparsity pattern with the nonzero values filled in) is returned.
 
"""
function makematrix!(self::TPATT) where {TPATT<:SysmatAssemblerSparsePattwLookup}
    return self._pattern
end
