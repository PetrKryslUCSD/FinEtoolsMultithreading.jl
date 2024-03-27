var documenterSearchIndex = {"docs":
[{"location":"guide/guide.html#Guide","page":"How to guide","title":"Guide","text":"","category":"section"},{"location":"guide/guide.html","page":"How to guide","title":"How to guide","text":"The FinEtools  package is used here to solve a variety of finite element problems. This package can provide an overlay to parallelize sparse matrix assembly.","category":"page"},{"location":"guide/guide.html","page":"How to guide","title":"How to guide","text":"There is one high-level function that can be used to parallelize the assembly of any sparse matrix – acoustic mass or stiffness,  conductivity matrix, stiffness or mass matrix, etc.","category":"page"},{"location":"guide/guide.html","page":"How to guide","title":"How to guide","text":"This bit of code will assemble the diffusion bilinear form sparse matrix stored in the CSR format:","category":"page"},{"location":"guide/guide.html","page":"How to guide","title":"How to guide","text":"function createsubdomain(fessubset)\n    FEMMBase(IntegDomain(fessubset, GaussRule(3, 2)))\nend\n\nfunction matrixcomputation!(femm, assblr)\n    bilform_diffusion(femm, assblr, geom, psi, DataCache(Matrix(1.0 * LinearAlgebra.I(3))))\nend\n\nK1 = parallel_make_matrix(fes, psi, createsubdomain, matrixcomputation!;\n    ntasks=ntasks, kind=:CSR)","category":"page"},{"location":"guide/guide.html","page":"How to guide","title":"How to guide","text":"The user can ask for any number of tasks to be used (even though it would be best to match it to the number of available threads).","category":"page"},{"location":"index.html#FinEtoolsMultithreading-Documentation","page":"Home","title":"FinEtoolsMultithreading Documentation","text":"","category":"section"},{"location":"index.html","page":"Home","title":"Home","text":"","category":"page"},{"location":"index.html#Conceptual-guide","page":"Home","title":"Conceptual guide","text":"","category":"section"},{"location":"index.html","page":"Home","title":"Home","text":"The construction of the toolkit is described: the composition of modules, the basic data structures, the methodology of computing quantities required in the finite element methodology, and more.","category":"page"},{"location":"index.html","page":"Home","title":"Home","text":"Pages = [\n    \"guide/guide.md\",\n]\nDepth = 1","category":"page"},{"location":"index.html#Manual","page":"Home","title":"Manual","text":"","category":"section"},{"location":"index.html","page":"Home","title":"Home","text":"The description of the types and the functions, organized by module and/or other logical principle.","category":"page"},{"location":"index.html","page":"Home","title":"Home","text":"Pages = [\n    \"man/man.md\",\n]\nDepth = 2","category":"page"},{"location":"index.html#Index","page":"Home","title":"Index","text":"","category":"section"},{"location":"index.html","page":"Home","title":"Home","text":"","category":"page"},{"location":"man/man.html#Manual","page":"Manual","title":"Manual","text":"","category":"section"},{"location":"man/man.html","page":"Manual","title":"Manual","text":"CurrentModule = FinEtoolsMultithreading","category":"page"},{"location":"man/man.html#High-level-API","page":"Manual","title":"High-level API","text":"","category":"section"},{"location":"man/man.html","page":"Manual","title":"Manual","text":"parallel_make_matrix","category":"page"},{"location":"man/man.html#FinEtoolsMultithreading.parallel_make_matrix","page":"Manual","title":"FinEtoolsMultithreading.parallel_make_matrix","text":"parallel_make_matrix(\n    fes,\n    u,\n    createsubd,\n    matrixupdt!;\n    ntasks = Threads.nthreads(),\n    kind = :CSC,\n)\n\nAssemble a sparse matrix.\n\nEither a :CSC matrix or a :CSR matrix is created. We shall refer to this matrix as a CSX matrix. The process is:\n\nConstruct the incidence relation node-to-neighbors.\nMake the sparse pattern and create a sparse CSX matrix with all values zero.\nConstruct the incidence relation element-to-neighbors.\nCompute element coloring.\nSet up domain decomposition.\nCompute and assemble the matrix entries.\n\n\n\n\n\nparallel_make_matrix(\n    fes,\n    dofnums,\n    ndofs,\n    FT,\n    n2e,\n    createsubd,\n    matrixupdt!,\n    ntasks,\n    kind,\n)\n\nAssemble a sparse matrix.\n\nConstruct the incidence relation node-to-neighbors.\nMake the sparse pattern and create a sparse CSX matrix with all values zero.\nConstruct the incidence relation element-to-neighbors.\nCompute element coloring.\nSet up domain decomposition.\nCompute and assemble the matrix entries.\n\n\n\n\n\n","category":"function"}]
}
