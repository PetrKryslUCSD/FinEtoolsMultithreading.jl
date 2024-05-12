using ECLGraphColor: PECLgraph, make_graph, add_nlist_all_row, add_nindex
using ECLGraphColor: get_color, run_graph_coloring, free_graph, print_graph, write_graph

"""
    element_coloring(fes, e2e::E2EM, ntasks::IT1) where {E2EM<:FElemToNeighborsMap{IT} where {IT},IT1<:Integer}
    
Determine element coloring such that no elements of the same color share a node.

# Arguments
- `fes`: The finite element set.
- `e2e`: The element-to-element map.

# Returns
Vector of element colors, vector of unique colors.
"""
function element_coloring(fes, e2e::E2EM, ntasks::IT1) where {E2EM<:FElemToNeighborsMap{IT} where {IT},IT1<:Integer}
    element_colors = fill(zero(Int16), count(fes))
    if !((Base.Sys.islinux()) || (Base.Sys.isapple())) 
        return element_coloring(fes, e2e)
    end
    map = e2e.map
    g = make_graph(length(map), sum([length(c) for c in map]))
    idx = 1
    for i in eachindex(map)
        add_nindex(g, i, idx)
        idx += length(map[i]) 
    end
    add_nindex(g, length(map)+1, idx)
    Threads.@threads for i in eachindex(map)
        neighbors = map[i]
        add_nlist_all_row(g, i, length(neighbors), neighbors); 
    end
    # print_graph(g)
    run_graph_coloring(g, ntasks, 0, 0)
    Threads.@threads for i in 1:length(map)
        element_colors[i] = get_color(g, i) 
    end
    # write_graph(g, "testgraph.egr")
    free_graph(g)
    return element_colors, sort(unique(element_colors))
end
