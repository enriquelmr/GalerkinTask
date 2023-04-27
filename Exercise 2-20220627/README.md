# Exercise 2.3 for "Numerical Methods for PDEs - Galerkin Methods", Summer Term 2022

First, you need to create the meshes as described in `disc_with_hole.geo`.
The file `disc_with_hole_standard_kernel.geo` sets up the same geometry using
the older standard kernel of GMSH instead of the newer OpenCASCADE kernel
(where the syntax of creating circles is different).

Then, you can reproduce the results by executing
```julia
include("exercise03_3.jl")
save_visualization("disc_with_hole_2.msh")
```
This will save files `model_0.vtu`, `model_1.vtu`, `model_2.vtu` where
`model_d.vtu` contains `d`-dimensional information. In particular,
`model_1.vtu` contains 1D information (of all surfaces of the triangles).
You can visualize the files for example with Paraview. Alternatively, you
can of course create the meshes and visualize them directly in GMSH.

Please see comments in the source files for further remarks.

The code was developed for Julia v1.7.3.
