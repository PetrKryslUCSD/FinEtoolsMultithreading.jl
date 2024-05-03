"""
    element_coloring(fes, e2e::E2EM, ellist::Vector{IT} = Int[]) where {E2EM<:FElemToNeighborsMap,IT<:Integer}
    
Determine element coloring such that no elements of the same color share a node.

# Arguments
- `fes`: The finite element set.
- `e2e`: The element-to-element map.
- `ellist`: Optional list of element serial numbers to include (default is all
  elements need to be considered).

# Returns
Vector of element colors, vector of unique colors, and vector of counts of each
color.

"""
function element_coloring(fes, e2e::E2EM, ellist::Vector{IT} = Int[]) where {E2EM<:FElemToNeighborsMap,IT<:Integer}
    element_colors = fill(zero(Int16), count(fes))
    unique_colors = eltype(element_colors)[1]
    color_counts = IT[0]
    if isempty(ellist) 
        ellistit = 1:count(fes)
    else
        ellistit = ellist
    end
    return element_coloring!(element_colors, unique_colors, color_counts, e2e.map, ellistit)
end

function __find_of_min_count(color_used, first, color_counts)
    c = 0
    mincount = typemax(eltype(color_counts))
    @inbounds for k in first:length(color_used)
        if color_used[k] == 0 && mincount > color_counts[k]
            mincount = color_counts[k]
            c = k
        end
    end
    return c
end

function __find_first_avail(color_used)
    @inbounds for k in eachindex(color_used)
        if color_used[k] == 0 
            return k
        end
    end
    return 0
end

function element_coloring!(element_colors, unique_colors, color_counts, e2emap, ellist_iterator) 
    color_used = fill(zero(eltype(element_colors)), length(unique_colors))
    for e in ellist_iterator
        if element_colors[e] == 0
            color_used .= 0; nchanges = 0
            for oe in e2emap[e]
                c = element_colors[oe]
                if c > 0
                    color_used[c] += 1
                    nchanges += 1
                end
            end
            if nchanges == 0
                c = __find_of_min_count(color_used, 1, color_counts)
                element_colors[e] = unique_colors[c]
                color_counts[c] += 1
            else
                first = __find_first_avail(color_used)
                if first == 0
                    added = maximum(unique_colors) + 1
                    push!(unique_colors, added)
                    push!(color_counts, 1)
                    push!(color_used, 0)
                    element_colors[e] = added
                else
                    c = __find_of_min_count(color_used, first, color_counts)
                    element_colors[e] = unique_colors[c]
                    color_counts[c] += 1
                end
            end
        end
    end
    return element_colors, unique_colors, color_counts
end