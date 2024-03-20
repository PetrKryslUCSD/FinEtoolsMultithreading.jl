
function _zeros_via_calloc(::Type{T}, dims::Integer...) where {T}
    ptr = Ptr{T}(Libc.calloc(prod(dims), sizeof(T)))
    return unsafe_wrap(Array{T}, ptr, dims; own = true)
end
