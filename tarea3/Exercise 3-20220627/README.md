# Exercise 3.3 for "Numerical Methods for PDEs - Galerkin Methods", Summer Term 2022

First, you need to create the meshes as described in `disc_with_hole.geo`
from Exercise 2.3 (the previous exercise sheet).

Then, you can reproduce the results by executing
```julia
include("exercise03_3.jl")
exercise03_3()
```
This will compute numerical solutions on each mesh in this
directory and save them to files `disc_with_hole_*_solution.vtu`.
These files can be opened with Paraview for manual inspection.
Moreover, a table of discrete L2 and H1 errors will be printed
to the screen.

By default, this uses piecewise linear finite elements. For quadratic
Lagrange finite elements, run
```julia
include("exercise03_3.jl")
exercise03_3(order=2)
```

Please see comments in the source files for further remarks.

The code was developed for Julia v1.7.3.
