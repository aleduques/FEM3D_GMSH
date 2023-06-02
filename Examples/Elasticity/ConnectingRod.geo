Merge "ConnectingRod.step";
Physical Volume("Connecting Rod", 31) = {1};
Physical Surface("Fixed Support", 32) = {9};
Physical Surface("Loaded Surface", 34) = {10};
Physical Surface("Free Surface", 35) = {12, 4, 3, 5, 6, 7, 8, 2, 1,11};


\\ Change element Size
Mesh.MeshSizeMax = 0.001;

\\ Change Element order
Mesh.ElementOrder = 2;
Mesh 3;