# FinEtoolsMultithreading.jl examples

The assumption is that the user is working in the `FinEtoolsMultithreading/examples` folder.
Furthermore, the use of a shell is presumed (for instance `bash`).

## Linear elasticity test case

To run the parallel and serial assembly for a mesh with 40 mesh edges along one
of the dimensions, type the following into the shell prompt:

The parallel assembly of the sparse matrix on three computing threads
```
julia -t 3 ./testlindef/parsim.jl 40 
```

The serial assembly of the sparse matrix 
```
julia ./testlindef/sersim.jl 40 
```