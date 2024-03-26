function element_coloring(fes, e2e::E2EM, ellist::Vector{IT} = Int[]) where {E2EM<:FEElemToNeighborsMap,IT<:Integer}
    element_colors = fill(zero(eltype(fes.conn[1])), count(fes))
    unique_colors = eltype(element_colors)[1]
    color_counts = eltype(element_colors)[0]
    if isempty(ellist) 
        ellistit = 1:count(fes)
    else
        ellistit = ellist
    end
    return element_coloring!(element_colors, unique_colors, color_counts, e2e.map, ellistit)
end

function __find_of_min_count(color_used, color_counts)
    c = 0
    mincount = typemax(eltype(color_counts))
    @inbounds for k in eachindex(color_used)
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
                c = __find_of_min_count(color_used, color_counts)
                element_colors[e] = unique_colors[c]
                color_counts[c] += 1
            else
                if __find_first_avail(color_used) == 0
                    added = maximum(unique_colors) + 1
                    push!(unique_colors, added)
                    push!(color_counts, 0)
                    push!(color_used, 0)
                    element_colors[e] = added
                    color_counts[c] += 1
                else
                    c = __find_of_min_count(color_used, color_counts)
                    element_colors[e] = unique_colors[c]
                    color_counts[c] += 1
                end
            end
        end
    end
    return element_colors, unique_colors, color_counts
end