// Gmsh project created on Fri Jun 24 11:25:32 2022
SetFactory("OpenCASCADE");
//+
Circle(1) = {0, 0, 0, 0.5, 0, 2*Pi};
//+
Circle(2) = {0, 0, 0, 2, 0, 2*Pi};
//+
Physical Curve("inner", 3) = {1};
//+
Physical Curve("outer", 4) = {2};
//+
Curve Loop(1) = {1};
//+
Curve Loop(2) = {2};
//+
Plane Surface(1) = {1, 2};
//+
Physical Surface("mesh", 5) = {1};
