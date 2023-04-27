// Gmsh project created on Thu Jun 23 12:14:01 2022
SetFactory("OpenCASCADE");
//+
Rectangle(1) = {0, 0, 0, 1, 1, 0};
//+
Physical Curve("wall", 5) = {3, 2, 1, 4};
//+
Physical Curve("wall", 5) += {3};
//+
Physical Curve("wall", 5) += {4, 3, 2, 1};
//+
Physical Surface("malla", 6) = {1};
//+
Show "*";
