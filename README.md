[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![Build status](https://github.com/PetrKryslUCSD/FinEtoolsMultithreading.jl/workflows/CI/badge.svg)](https://github.com/PetrKryslUCSD/FinEtoolsMultithreading.jl/actions)
[![Code Coverage](https://codecov.io/gh/PetrKryslUCSD/FinEtoolsMultithreading.jl/branch/master/graph/badge.svg)](https://app.codecov.io/gh/PetrKryslUCSD/FinEtoolsMultithreading.jl)
[![Latest documentation](https://img.shields.io/badge/docs-latest-blue.svg)](https://petrkryslucsd.github.io/FinEtoolsMultithreading.jl/latest)



# FinEtoolsMultithreading.jl

This package is an overlay  for FinEtools-based application packages. It can be
used to parallelize certain finite element operations with multithreading.

Currently, sparse matrix operators can be assembled in parallel on multiple threads.
