"""
    parallel_element_coloring(fes, e2e::E2EM, ellist::Vector{IT} = Int[]) where {E2EM<:FElemToNeighborsMap,IT<:Integer}
    
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
function parallel_element_coloring(fes, e2e::E2EM,
    ellist::Vector{IT}=Int[]) where {E2EM<:FElemToNeighborsMap,IT<:Integer}
    element_colors = fill(zero(Int16), count(fes))
    if isempty(ellist)
        ellistit = 1:count(fes)
    else
        ellistit = ellist
    end
    return parallel_element_coloring!(element_colors, e2e.map, ellistit)
end

function _extrema(eneighbors, element_colors)
    mine = typemax(eltype(eneighbors))
    maxe = zero(eltype(eneighbors))
    for oe in eneighbors
        element_colors[oe] == 0 && (mine = min(mine, oe))
        element_colors[oe] == 0 && (maxe = max(maxe, oe))
    end
    return mine, maxe
end

function _minimum(eneighbors, element_colors)
    mine = typemax(eltype(eneighbors))
    @inbounds for oe in eneighbors
        element_colors[oe] == 0 && (mine = min(mine, oe))
    end
    return mine
end

function _avail_color(eneighbors, element_colors, current_color)
    color_used = fill(false, current_color)
    @inbounds for oe in eneighbors
        c = element_colors[oe]
        if c != 0
            color_used[c] = true
        end
    end
    @inbounds for c in eachindex(color_used)
        if !color_used[c]
            return c
        end
    end
    return 0
end

# function parallel_element_coloring!(element_colors, e2emap, ellist_iterator)
#     total_elements = length(ellist_iterator)
#     new_element_colors = deepcopy(element_colors)
#     colored_elements = 0
#     current_color = 1
#     done = false
#     while !done
#         increase_min_color = false
#         increase_max_color = false
#         for e in ellist_iterator
#             if element_colors[e] == 0 # not colored yet
#                 mine, maxe = _extrema(e2emap[e], element_colors)
#                 if (e == mine)
#                     c = _avail_color(e2emap[e], element_colors, current_color-1)
#                     if c == 0
#                         c = current_color
#                         increase_min_color = true
#                     end 
#                     @assert !(c in element_colors[e2emap[e]])
#                     @assert !(c in new_element_colors[e2emap[e]])
#                     new_element_colors[e] = c
#                     colored_elements += 1
#                     # @show  e, new_element_colors[e], element_colors[e2emap[e]], new_element_colors[e2emap[e]]
#                 elseif (e == maxe)
#                     @show c = _avail_color(e2emap[e], element_colors, current_color-1)
#                     if c == 0
#                         c = current_color + 1
#                         increase_max_color = true
#                     end 
#                     @assert !(c in element_colors[e2emap[e]])
#                     if (c in new_element_colors[e2emap[e]])
#                         @show c, new_element_colors[e2emap[e]]
#                     end
#                     @assert !(c in new_element_colors[e2emap[e]])
#                     new_element_colors[e] = c
#                     colored_elements += 1
#                     # @show  e, new_element_colors[e], element_colors[e2emap[e]], new_element_colors[e2emap[e]]
#                 end
                   
#             end
#         end
#         increase_min_color && (current_color += 1)
#         increase_max_color && (current_color += 1)
#         @show current_color
#         @show colored_elements, total_elements
#         element_colors .= new_element_colors 
#         done = colored_elements == total_elements
#     end
#     return element_colors, collect(1:current_color-1)
# end

# function parallel_element_coloring!(element_colors, e2emap, ellist_iterator)
#     total_elements = length(ellist_iterator)
#     new_element_colors = deepcopy(element_colors)
#     colored_elements = 0
#     current_color = 1
#     done = false
#     while !done
#         increase_current_color = false
#         color_used = fill(false, current_color-1)
#         new_ellist_iterator = eltype(ellist_iterator)[]
#         sizehint!(new_ellist_iterator, length(ellist_iterator))
#         for e in ellist_iterator
#             if element_colors[e] == 0 # not colored yet
#                 mine = _minimum(e2emap[e], element_colors)
#                 if (e == mine)
#                     c = _avail_color(e2emap[e], element_colors, color_used)
#                     if c == 0
#                         c = current_color
#                         increase_current_color = true
#                     end 
#                     new_element_colors[e] = c
#                     colored_elements += 1
#                 else
#                     push!(new_ellist_iterator, e)
#                 end
#             end
#         end
#         increase_current_color && (current_color += 1)
#         # @show current_color
#         # @show colored_elements, total_elements
#         element_colors .= new_element_colors 
#         ellist_iterator = new_ellist_iterator
#         done = colored_elements == total_elements
#     end
#     return element_colors, collect(1:current_color-1)
# end

# function parallel_element_coloring!(element_colors, e2emap, ellist_iterator)
#     total_elements = length(ellist_iterator)
#     new_element_colors = deepcopy(element_colors)
#     colored_elements = 0
#     current_color = 1
#     done = false
#     while !done
#         increase_current_color = false
#         color_used = fill(false, current_color-1)
#         Threads.@threads for e in ellist_iterator
#             if element_colors[e] == 0 # not colored yet
#                 mine = _minimum(e2emap[e], element_colors)
#                 if (e == mine)
#                     c = _avail_color(e2emap[e], element_colors, color_used)
#                     if c == 0
#                         c = current_color
#                         increase_current_color = true
#                     end 
#                     new_element_colors[e] = c
#                     colored_elements += 1
#                 end
#             end
#         end
#         increase_current_color && (current_color += 1)
#         element_colors .= new_element_colors 
#         done = colored_elements == total_elements
#     end
#     return element_colors, collect(1:current_color-1)
# end

# function parallel_element_coloring!(element_colors, e2emap, ellist_iterator)
#     total_elements = length(ellist_iterator)
#     new_element_colors = deepcopy(element_colors)
#     colored_elements = 0
#     current_color = 1
#     done = false
#     while !done
#         increase_current_color = false
#         @inbounds for e in ellist_iterator
#             if element_colors[e] == 0 # not colored yet
#                 mine = _minimum(e2emap[e], element_colors)
#                 if (e == mine)
#                     c = _avail_color(e2emap[e], element_colors, current_color-1)
#                     if c == 0
#                         c = current_color
#                         increase_current_color = true
#                     end 
#                     new_element_colors[e] = c
#                     colored_elements += 1
#                 end
#             end
#         end
#         increase_current_color && (current_color += 1)
#         element_colors .= new_element_colors 
#         done = colored_elements == total_elements
#     end
#     return element_colors, collect(1:current_color-1)
# end

function _maybe_color_element!(new_element_colors, e, eneighbors, element_colors, current_color)
    mine = _minimum(eneighbors, element_colors)
    if (e == mine)
        increase_current_color = false
        c = _avail_color(eneighbors, element_colors, current_color - 1)
        if c == 0
            c = current_color
            increase_current_color = true
        end
        new_element_colors[e] = c
        return increase_current_color, 1
    else
        return false, 0
    end
end

function parallel_element_coloring!(element_colors, e2emap, ellist_iterator)
    ntasks = Threads.nthreads()
    chks = chunks(ellist_iterator, ntasks)
    new_element_colors = deepcopy(element_colors)
    increase_current_color = fill(false, ntasks)
    colored_elements = fill(zero(eltype(ellist_iterator)), ntasks)
    current_color = 1
    remaining = length(ellist_iterator)
    while remaining > 0
        increase_current_color .= false
        colored_elements .= 0
        Threads.@sync begin
            Threads.@spawn for ch in chks
                t = ch[2]
                for i in ch[1]
                    e = ellist_iterator[i]
                    if element_colors[e] == 0 # not colored yet
                        inc, col = _maybe_color_element!(new_element_colors, e, e2emap[e], element_colors, current_color)
                        increase_current_color[t] = increase_current_color[t] || inc
                        colored_elements[t] += col
                    end
                end
            end
        end
        any(increase_current_color) && (current_color += 1)
        remaining -= sum(colored_elements)
        element_colors .= new_element_colors 
    end
    return element_colors, collect(1:current_color-1)
end
