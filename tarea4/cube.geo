// Gmsh project created on Mon Jun 27 16:44:58 2022
SetFactory("OpenCASCADE");
//+
Box(1) = {0, 0, 0, 1, 1, 1};
//+
Physical Curve("wall", 13) = {11, 12, 3, 4, 2, 7, 8, 6, 9, 10, 1, 5};
//+
Physical Volume("mesh3d", 14) = {1};
