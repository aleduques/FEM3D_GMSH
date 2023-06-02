# FEM3D_GMSH

This repository contains a vectorized Matlab/Octave package for solving three-dimensional finite element problems with meshes obtained from GMSH. Currently, it supports tetrahedral elements up to degree 4.


### Installation instructions
Follow these steps to install the package:

1. Navigate to the "code" section of the repository and click on "Download ZIP" to obtain all the necessary files. Save the downloaded files in a folder.

2. To use the package, open the command line and type "__addpath folderDir__," where "folderDir" is the directory path to the folder containing the package.
### Documentation

This code was developed for a final project in Industrial Engineering. For more __detailed information__ on the class and basic ideas on the Finite Element Method, please refer to the bachelor's thesis: 
[FEM3D_simplicialGMSH.pdf](https://github.com/aleduques/FEM3D_GMSH/blob/main/TFG_ALEJANDRO_DUQUE_SALAZAR.pdf)


### Current Folders:
1.  FEM3Dclass:  This folder contains the class with various available functions:
  - BasisNj: Lagrange Basis in the tetrahedral reference element and their partial derivatives. Lagrange basis in the            triangular element.
  - evalFEM3DUh: Compute Finite Element Solution at non-nodal points
  - femAdvectionMatrix, femMassMatrix, femRobin, femSourceTerm, femStiffnessMatrix, femStressMatrix:  Assembly system of        equations.
  - quadRule2D and quadRule3D: Quadrature formulas valid for the tetrahedral and  triangular reference element with nodes        given in barycentric coordinates

  - local3DMatrices and localRobinMass: Exact local matrices.


2.  Examples:  This folder contains two separate subfolders, one for Linear Elasticity and the other for the Heat Equation. Each subfolder includes a .m script (Matlab and Octave) and a mesh file created in GMSH. Additionally, the corresponding .geo file used to create the meshes is provided. In this case, the meshes were created with the assistance of FreeCAD, and the geometries were stored in files with the .step extension.


