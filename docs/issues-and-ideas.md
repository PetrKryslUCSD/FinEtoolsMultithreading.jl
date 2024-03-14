
Issues and ideas:


-- Documenter:
using FinEtoolsMultithreading
using DocumenterTools
Travis.genkeys(user="PetrKryslUCSD", repo="https://github.com/PetrKryslUCSD/FinEtoolsMultithreading.jl")

using Pkg; Pkg.add("DocumenterTools");                                 
using DocumenterTools                                                  
DocumenterTools.genkeys(user="PetrKryslUCSD", repo="git@github.com:PetrKryslUCSD/FinEtoolsMultithreading.jl.git")                                                  
using Pkg; Pkg.rm("DocumenterTools");  
