// Create a mesh of a disc with a hole in the center
//
// You can adapt the mesh parameter globally by the command line argument
//   `-clscale FACTOR`
// of GMSH. For example, to create a series of refined meshes, you can run
//   gmsh disc_with_hole.geo -2 -clscale 1 -o disc_with_hole_1.msh
//   gmsh disc_with_hole.geo -2 -clscale 0.5 -o disc_with_hole_2.msh
//   gmsh disc_with_hole.geo -2 -clscale 0.25 -o disc_with_hole_3.msh
//   gmsh disc_with_hole.geo -2 -clscale 0.125 -o disc_with_hole_4.msh
//   gmsh disc_with_hole.geo -2 -clscale 0.0625 -o disc_with_hole_5.msh

// Use the OpenCASCADE kernel. This changes the syntax of some commands below.
SetFactory("OpenCASCADE");

// Boundary curves
// Circle(TAG) = {center_x, center_y, center_z, radius, angle_min, angle_max};
// First, we create "elementary curves"
// for the inner circle
Circle(1) = {0.0, 0.0, 0.0, 0.5, 0.0, 2*Pi};
// and the outer circle.
Circle(2) = {0.0, 0.0, 0.0, 2.0, 0.0, 2*Pi};
// Next, we create "curve loops" for GMSH.
Curve Loop(3) = {1};
Curve Loop(4) = {2};

// Create the surface with hole
// Syntax: Outer curve loop, inner curve loops for holes
Plane Surface(5) = {4, 3};

// Now, we create physical groups with names.
// NOTE: The `Circle` commands above create an additional point numbered
//       automatically. Here, these numbers are the same as the ones of the
//       corresponding circles. These points are by default not added to the
//       physical boundaries when adding the curves. Thus, we need to add them
//       manually.
Physical Curve("inner_boundary") = {1};
Physical Point("inner_boundary") = {1};
Physical Curve("outer_boundary") = {2};
Physical Point("outer_boundary") = {2};
Physical Surface("domain") = {5};
