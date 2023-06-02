Merge "cylinderFins.step";


// Define Physical Volume
Physical Volume("Finned Cylinder", 1) = {1};


// Define Physical Surfaces
interior = {38};
exterior[] = Surface{:};
exterior[] -= interior[];

Physical Surface("Interior") = {interior[]};
Physical Surface("Exterior") = {exterior[]};

// Mesh
Mesh.MeshSizeMax = 0.1;
Mesh 3;


