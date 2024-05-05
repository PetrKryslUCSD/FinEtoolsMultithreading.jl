module mparallelcoloring_2
using FinEtools
using FinEtoolsMultithreading.Exports
using FinEtoolsMultithreading: parallel_element_coloring
using LinearAlgebra
using Test

function test_coloring(coloring, n2e)
    element_colors, unique_colors = coloring
    findall(c -> c <= 0, element_colors) === nothing
    findall(c -> c > maximum(unique_colors), element_colors) === nothing
    # @show sort(unique(element_colors)), sort(unique_colors)
    @assert norm(sort(unique(element_colors)) - sort(unique_colors)) == 0
    for k in eachindex(n2e.map)
        nc = element_colors[n2e.map[k]]
        @assert length(nc) == length(unique(nc))
        @assert norm(sort(nc) - sort(unique(nc))) == 0
    end
    true
end

function test()
    W = 11.0
    L = 12.0
    t = 10.0
    nl, nt, nw = 72, 83, 24
    nl, nt, nw = 2, 3, 4
    # nl, nt, nw = 7, 13, 9
    # nl, nt, nw = 17, 13, 19
    # nl, nt, nw = 27, 23, 29
    ntasks = 2

    fens, fes = T4block(Float64.((nl, nw, nt))..., nl, nw, nt)
    
    n2e = FENodeToFEMap(fes, count(fens))
    n2n = FENodeToNeighborsMap(n2e, fes.conn)
    e2e = FElemToNeighborsMap(n2e, fes.conn)
    
    coloring = element_coloring(fes, n2e, collect(1:count(fes)))
    @time coloring = element_coloring(fes, n2e, collect(1:count(fes)))
    @show unique(coloring[1])
    test_coloring(coloring, n2e)

    coloring = parallel_element_coloring(fes, e2e, collect(1:count(fes)))
    @time coloring = parallel_element_coloring(fes, e2e, collect(1:count(fes)))
    @show unique(coloring[1])
 

    element_colors, unique_colors = coloring
    @show unique_colors = sort(unique_colors)
    @show partitions = unique_colors
    for j in 1:length(partitions)
        sfes = subset(fes, findall(v -> v == partitions[j], element_colors))
        @show count(sfes)
        vtkexportmesh("mesh_test_coloring_1-c=$(partitions[j]).vtk", fens, sfes)
    end
    
   test_coloring(coloring, n2e)

    true
end
test()
nothing
end