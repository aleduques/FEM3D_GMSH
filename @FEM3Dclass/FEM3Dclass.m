% Some explorations
%
% To initialize the class run
%
% T = FEM3Dclass(filename)
%
% where filename is a mesh from GMSH (4.1 version, ascii format)
%
% The class has the following properties
%
% mesh : mesh file
% Nj3D : the 3D local basis (unit tetrahedra)
% Nj2D : the 2D local basis (unit triangle)
%
% 04 01 2022 by vD-aD
%
% Last modification, the constructor is merely a prototype

classdef FEM3Dclass < handle
    properties
        mesh     =[];   % Mesh structure
        Nj3D     =[];   % local 3D basis
        Nj2D     =[];   % local 2D basis
        gradNj3D =[];
        ComplementaryInformation = [];
    end

    methods (Static)


        [T,S] = importGMSH3D(var)

        % Quadrature/Cubature rules
        quadRuleTtrh = quadRule3D
        quadRuletr   = quadRule2D

        % Local matrices
        [M,Sxx,Sxy,Sxz,Syy,Syz,Szz,Ax,Ay,Az] = local3DMatrices(deg)
        R = localRobinMass(deg)

        %  Lagrange basis/Gradients
        [Nj2D,Nj3D, gradNj3D] = basisNj(deg)

    end

    methods
        % constructor
        function obj = FEM3Dclass(var)

            if ischar(var) && exist(var, 'file') == 2 %Mesh file
                [T, ComplementaryInformation] = FEM3Dclass.importGMSH3D(var);
                obj.mesh = T;
                obj.ComplementaryInformation = ComplementaryInformation;
            elseif isstruct(var) % Loaded mesh variable
                obj.mesh = var;
            else                   
            error('The file chosen does not exist')        
            end
            [Nj2D,Nj3D,gradNj3D] = FEM3Dclass.basisNj(obj.mesh.deg);
            obj.Nj2D = Nj2D;
            obj.Nj3D = Nj3D;
            obj.gradNj3D = gradNj3D;
        end

        % Assembly
        M = femMassMatrix(obj, c)

        S = femStressMatrix(obj)

        S = femStiffnessMatrix(obj, K)

        b = femSourceTerm(obj, f)

        [t,MR] = femRobin(obj, robin, gn, varargin)

        [A] =  femAdvectionMatrix(obj, v)
        
        
        % Finite element evaluation
        [val,Meval,barPt,indTtrh,indPtError] = evalFEM3DUh(uh,obj,pt,varargin)


    end


end


