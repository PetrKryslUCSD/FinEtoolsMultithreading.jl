using ECLGraphColor: PECLgraph, make_graph, add_nlist_all_row, add_nindex
using ECLGraphColor: get_color, run_graph_coloring, free_graph, print_graph, write_graph

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
    ellist::Vector{IT}=Int[], ntasks = Threads.nthreads()) where {E2EM<:FElemToNeighborsMap,IT<:Integer}
    element_colors = fill(zero(Int16), count(fes))
    if (!Base.Sys.islinux())
        return element_coloring(fes, e2e, ellist)
    end
    if isempty(ellist)
        ellistit = 1:count(fes)
    else
        ellistit = ellist
    end
    map = e2e.map
@time begin
    g = make_graph(length(map), sum([length(c) for c in map]))
    idx = 1
    for i in eachindex(map)
        add_nindex(g, i, idx)
        idx += length(map[i]) 
    end
    add_nindex(g, length(map)+1, idx)
end
@time begin
    Threads.@threads for i in eachindex(map)
        neighbors = map[i]
        add_nlist_all_row(g, i, length(neighbors), neighbors); 
    end
end
    # print_graph(g)
    run_graph_coloring(g, ntasks, 0, 1)
@time begin
    Threads.@threads for i in 1:length(map)
        element_colors[i] = get_color(g, i) 
    end
end
    # write_graph(g, "testgraph.egr")
    free_graph(g)
    return element_colors, sort(unique(element_colors))
end

# function parallel_element_coloring(fes, e2e::E2EM,
#     ellist::Vector{IT}=Int[]) where {E2EM<:FElemToNeighborsMap,IT<:Integer}
#     element_colors = fill(zero(Int16), count(fes))
#     if isempty(ellist)
#         ellistit = 1:count(fes)
#     else
#         ellistit = ellist
#     end
#     return parallel_element_coloring!(element_colors, e2e.map, ellistit)
# end

# function _extrema(eneighbors, element_colors)
#     mine = typemax(eltype(eneighbors))
#     maxe = zero(eltype(eneighbors))
#     for oe in eneighbors
#         element_colors[oe] == 0 && (mine = min(mine, oe))
#         element_colors[oe] == 0 && (maxe = max(maxe, oe))
#     end
#     return mine, maxe
# end

# function _minimum(eneighbors, element_colors)
#     mine = typemax(eltype(eneighbors))
#     @inbounds for oe in eneighbors
#         element_colors[oe] == 0 && (mine = min(mine, oe))
#     end
#     return mine
# end

# function _avail_color(eneighbors, element_colors, current_color)
#     color_used = fill(false, current_color)
#     @inbounds for oe in eneighbors
#         c = element_colors[oe]
#         if c != 0
#             color_used[c] = true
#         end
#     end
#     @inbounds for c in eachindex(color_used)
#         if !color_used[c]
#             return c
#         end
#     end
#     return 0
# end


# Working sequential  implementation
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

# function _maybe_color_element!(new_element_colors, e, eneighbors, element_colors, current_color)
#     mine = _minimum(eneighbors, element_colors)
#     if (e == mine)
#         increase_current_color = false
#         c = _avail_color(eneighbors, element_colors, current_color - 1)
#         if c == 0
#             c = current_color
#             increase_current_color = true
#         end
#         new_element_colors[e] = c
#         return increase_current_color, 1
#     else
#         return false, 0
#     end
# end

# # Working parallel implementation
# function parallel_element_coloring!(element_colors, e2emap, ellist_iterator)
#     ntasks = Threads.nthreads()
#     chks = chunks(ellist_iterator, ntasks)
#     new_element_colors = deepcopy(element_colors)
#     increase_current_color = fill(false, ntasks)
#     colored_elements = fill(zero(eltype(ellist_iterator)), ntasks)
#     current_color = 1
#     remaining = length(ellist_iterator)
#     while remaining > 0
#         increase_current_color .= false
#         colored_elements .= 0
#         Threads.@threads for ch in chks
#             t = ch[2]
#             for e in ellist_iterator[ch[1]]
#                 if element_colors[e] == 0 # not colored yet
#                     inc, col = _maybe_color_element!(new_element_colors, e, e2emap[e], element_colors, current_color)
#                     increase_current_color[t] = increase_current_color[t] || inc
#                     colored_elements[t] += col
#                 end
#             end
#         end
#         any(increase_current_color) && (current_color += 1)
#         remaining -= sum(colored_elements)
#         element_colors .= new_element_colors 
#     end
#     return element_colors, collect(1:current_color-1)
# end
