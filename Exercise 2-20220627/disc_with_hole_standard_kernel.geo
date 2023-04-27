// Create a mesh of a disc with a hole in the center
//
// You can adapt the mesh parameter globally by the command line argument
//   `-clscale FACTOR`
// of GMSH. For example, to create a series of refined meshes, you can run
//   cz
//   gmsh disc_with_hole_standard_kernel.geo -2 -clscale 0.5 -o disc_with_hole_2.msh
//   gmsh disc_with_hole_standard_kernel.geo -2 -clscale 0.25 -o disc_with_hole_3.msh
//   gmsh disc_with_hole_standard_kernel.geo -2 -clscale 0.125 -o disc_with_hole_4.msh
//   gmsh disc_with_hole_standard_kernel.geo -2 -clscale 0.0625 -o disc_with_hole_5.msh

// Boundary curves
// Circle(TAG) = {first_point_on_arc, center_of_circle, second_point_on_arc};
Point(1) = {0.0, 0.0, 0.0}; // center of the circles

Point(2) = {0.5, 0.0, 0.0};
Point(3) = {0.0, 0.5, 0.0};
Point(4) = {-0.5, 0.0, 0.0};
Point(5) = {0.0, -0.5, 0.0};
Circle(1) = {2, 1, 3};
Circle(2) = {3, 1, 4};
Circle(3) = {4, 1, 5};
Circle(4) = {5, 1, 2};

Point(6) = {2.0, 0.0, 0.0};
Point(7) = {0.0, 2.0, 0.0};
Point(8) = {-2.0, 0.0, 0.0};
Point(9) = {0.0, -2.0, 0.0};
Circle(5) = {6, 1, 7};
Circle(6) = {7, 1, 8};
Circle(7) = {8, 1, 9};
Circle(8) = {9, 1, 6};

// Next, we create "curve loops" for GMSH.
Curve Loop(1) = {1, 2, 3, 4};
Curve Loop(2) = {5, 6, 7, 8};

// Create the surface with hole
// Syntax: Outer curve loop, inner curve loops for holes
Plane Surface(1) = {2, 1};

// Now, we create physical groups with names.
Physical Curve("inner_boundary") = {1, 2, 3, 4};
Physical Point("inner_boundary") = {2, 3, 4, 5};
Physical Curve("outer_boundary") = {5, 6, 7, 8};
Physical Curve("outer_boundary") = {6, 7, 8, 9};
Physical Surface("domain") = {1};
