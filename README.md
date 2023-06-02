# FEM3D_GMSH

Vectorized Matlab/Octave package for solving three dimensional Finite element problems with meshes coming from GMSH. At the moment it works for tetrahedral elements up to degree 4.

### Installation instructions

1. Go to code, click on Download ZIP and store all information provided in a folder.

2. To use the class write in the command Line "addpath folderDir" where folderDir is the directory of the Folder containing the class.

### Documentation

This code was done for a final project in Industrial engineering. For more information check out 
[FEM3D_simplicialGMSH.pdf](https://github.com/aleduques/FEM3D_GMSH/blob/main/TFG_ALEJANDRO_DUQUE_SALAZAR.pdf)


### Current Folders:
Currently we have the folder @FEM3Dclass containing the class. The functions currently available are:
- BasisNj: Lagrange Basis in the tetrahedral reference element and their partial derivatives. Lagrange basis in the triangular element.
- evalFEM3DUh: Compute Finite Element Solution at non-nodal points
- femAdvectionMatrix, femMassMatrix, femRobin, femSourceTerm, femStiffnessMatrix, femStressMatrix: Assembly system of equations.
- quadRule2D and quadRule3D: Quadrature formulas for the tetrahedral and triangular reference element with nodes given in barycentric coordinates

- local3DMatrices and localRobinMass: Exact local matrices.

