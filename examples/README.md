# FinEtoolsMultithreading.jl examples

The assumption is that the user is working in the `FinEtoolsMultithreading/examples` folder.
Furthermore, the use of a shell is presumed (for instance `bash`).

Then, we could use the shell commands:

```
user:~$ mkdir try
user:~$ cd try
user:~/try$ git clone https://github.com/PetrKryslUCSD/FinEtoolsMultithreading.jl.git
user:~/try$ cd FinEtoolsMultithreading.jl/examples/
```

Then, we pick a testcase and proceed as described below. It is assumed that the executable of `julia` is on the path.

## Linear elasticity test case

The mesh with 50 mesh edges amounts to 500,000 elements. 

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


## Modal acoustic analysis test case

The mesh with number of refinements 6 amounts to 8.4 million elements. 

To run the parallel and serial assembly for a mesh refinement 4, type the following into the shell prompt:

The parallel assembly of the sparse matrix on three computing threads
```
julia -t 3 ./testacoust/parsim.jl 4
```

The serial assembly of the sparse matrix 
```
julia ./testacoust/sersim.jl 4
```
