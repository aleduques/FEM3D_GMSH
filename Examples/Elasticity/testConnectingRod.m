%%% Displacement analysis of a connecting rod.

% Mesh file
% T = FEM3Dclass('ConnectingRod_order2.msh'); 
T = FEM3Dclass('ConnectingRod.msh'); % Comment lines 136-139. Uncomment 140
%% Geometry Visualization


Loaded = T.ComplementaryInformation('Loaded Surface');
Loaded = T.mesh.domBd==Loaded;
LoadedSurface = T.mesh.trB(Loaded, 1:3);

Fixed = T.ComplementaryInformation('Fixed Support');
Fixed = T.mesh.domBd==Fixed;
FixedSurface = T.mesh.trB(Fixed, 1:3);

Free = T.ComplementaryInformation('Free Surface');
Free = ismember(T.mesh.domBd,Free);
FreeSurface = T.mesh.trB(Free, 1:3);

%
figure()
hold on
trimesh(LoadedSurface, T.mesh.coord(:,1),T.mesh.coord(:,2), ...
    T.mesh.coord(:,3),'FaceColor', 'red', 'EdgeColor','black')
trimesh(FixedSurface, T.mesh.coord(:,1),T.mesh.coord(:,2), ...
    T.mesh.coord(:,3),'FaceColor', 'blue', 'EdgeColor','black')
trimesh(FreeSurface, T.mesh.coord(:,1),T.mesh.coord(:,2), ...
        T.mesh.coord(:,3),'FaceColor', '#808080', 'EdgeColor','black')

axis('equal')
view(3)
%% Data

%%% Material Properties
% Young Modulus
E = 2e11; %Pa
% Poisson Coefficient
nu = 0.3;
%%% Lame Parameters
lambda = E*nu/(1+nu)/(1-2*nu);
mu = E/(1+nu)/2;
%% Load
q = 20;
q = q*1e6; %MPa



% Director cosines
LoadedUnitNormal = T.mesh.trBNormal(Loaded,:);
LoadedUnitNormal = LoadedUnitNormal./vecnorm(LoadedUnitNormal,2,2);
cosThetaX = LoadedUnitNormal(:,1);
cosThetaY = LoadedUnitNormal(:,2);
cosThetaZ = LoadedUnitNormal(:,3);

% Components of the load vector
qX =@(x,y,z) -q.*cosThetaX + x.*0;
qY =@(x,y,z) -q.*cosThetaY + x.*0;
qZ =@(x,y,z) -q.*cosThetaZ + x.*0;

%% Boundary Condition nodes
% Fixed Nodes
fixedNodes    = T.mesh.trB(Fixed,:); fixedNodes = unique(fixedNodes);
% Load
loadedNodes = T.mesh.trB(Loaded,:); loadedNodes = unique(loadedNodes);


nNodes = T.mesh.nNodes;
% Fixed support. All components of the displacement are zero
iD = fixedNodes; iD2 = [iD; iD+nNodes; iD+nNodes*2];
% Non Dirichlet nodes
inD = (1:length(T.mesh.coord)*3)'; inD = setdiff(inD, iD2);

%% System Construction

%%% Assembly

% Stiffness Matrix
disp('System Assembly')
A = {@(x,y,z) mu+2*lambda+x.*0, @(x,y,z) 0+x.*0,  @(x,y,z) 0+x.*0;
     @(x,y,z) 0+x.*0,           @(x,y,z) mu+x.*0, @(x,y,z) 0+x.*0;
     @(x,y,z) 0+x.*0,           @(x,y,z) 0+x.*0,  @(x,y,z) mu+x.*0};
S11 = femStiffnessMatrix(T,A);

A = {@(x,y,z) mu+x.*0, @(x,y,z) 0+x.*0,           @(x,y,z) 0+x.*0;
     @(x,y,z) 0+x.*0,  @(x,y,z) lambda+2*mu+x.*0, @(x,y,z) 0+x.*0;
     @(x,y,z) 0+x.*0,  @(x,y,z) 0+x.*0,           @(x,y,z) mu+x.*0};
S22 = femStiffnessMatrix(T,A);

A = {@(x,y,z) mu+x.*0, @(x,y,z) 0+x.*0,  @(x,y,z) 0+x.*0;
     @(x,y,z) 0+x.*0,  @(x,y,z) mu+x.*0, @(x,y,z) 0+x.*0;
     @(x,y,z) 0+x.*0,  @(x,y,z) 0+x.*0,  @(x,y,z) lambda+2*mu+x.*0};
S33 = femStiffnessMatrix(T,A);

A = {@(x,y,z) 0+x.*0,  @(x,y,z) lambda+x.*0, @(x,y,z) 0+x.*0;
     @(x,y,z) mu+x.*0, @(x,y,z) 0+x.*0,      @(x,y,z) 0+x.*0;
     @(x,y,z) 0+x.*0,  @(x,y,z) 0+x.*0,      @(x,y,z) 0+x.*0};
S12 = femStiffnessMatrix(T,A);

A = {@(x,y,z) 0+x.*0,  @(x,y,z) 0+x.*0, @(x,y,z) lambda+x.*0;
     @(x,y,z) 0+x.*0,  @(x,y,z) 0+x.*0, @(x,y,z) 0+x.*0;
     @(x,y,z) mu+x.*0, @(x,y,z) 0+x.*0, @(x,y,z) 0+x.*0};
S13 = femStiffnessMatrix(T,A);

A = {@(x,y,z) 0+x.*0, @(x,y,z) 0+x.*0,  @(x,y,z) 0+x.*0;
     @(x,y,z) 0+x.*0, @(x,y,z) 0+x.*0,  @(x,y,z) lambda+x.*0;
     @(x,y,z) 0+x.*0, @(x,y,z) mu+x.*0, @(x,y,z) 0+x.*0};
S23 = femStiffnessMatrix(T,A);

matrix = [S11  S12  S13
          S12' S22  S23
          S13' S23' S33];

% Neglect body forces
F = zeros(nNodes*3,1);

% Load
tX = femRobin(T, Loaded, qX);
tY = femRobin(T, Loaded, qY);
tZ = femRobin(T, Loaded, qZ);
t = [tX; tY; tZ];

disp('done')
% Solution initialization
u = zeros(nNodes*3,1);

% System Resolution
disp('Solving System')
F = F+t; F = F(inD);
F = F-matrix(inD, iD2)*u(iD2);
matrix = matrix(inD, inD);
u(inD) = matrix\F;

%%
uX = u(1:nNodes);
uY = u(nNodes+1:2*nNodes);
uZ = u(2*nNodes+1:3*nNodes);
%% PLOT SOLUTION

%%% Valid for P2 elements
% SurfaceSubTriangles = [1 4 6; 2 5 4; 3 6 5; 6 4 5];
% 
% 
% trB = [T.mesh.trB(:,SurfaceSubTriangles(1,:));
%                  T.mesh.trB(:,SurfaceSubTriangles(2,:));
%                  T.mesh.trB(:,SurfaceSubTriangles(4,:))];
%                  T.mesh.trB(:,SurfaceSubTriangles(3,:));

%%% Valid for P1 elements
trB = T.mesh.trB;
                 
fig = figure;
trisurf(trB, T.mesh.coord(:,1),T.mesh.coord(:,2),T.mesh.coord(:,3),...
            vecnorm([uX uY uZ],2,2), 'FaceColor','interp','EdgeColor','none');
colorbar
colormap turbo
xlabel('x')
ylabel('y')
zlabel('z')
axis('equal')
view(45,30)
title('Displacement module (m)')

fig = figure;
trisurf(trB, T.mesh.coord(:,1),T.mesh.coord(:,2),T.mesh.coord(:,3),...
            uX, 'FaceColor','interp','EdgeColor','none');
colorbar
colormap turbo
xlabel('x')
ylabel('y')
zlabel('z')
axis('equal')
view(45,30)
title('Displacement uX(m)')
%
fig = figure;
trisurf(trB, T.mesh.coord(:,1),T.mesh.coord(:,2),T.mesh.coord(:,3),...
            uY, 'FaceColor','inter p','EdgeColor','none');
colorbar
colormap turbo
xlabel('x')
ylabel('y')
zlabel('z')
axis('equal')
view(45,30)
title('displacement uY(m)')



fig = figure;
trisurf(trB, T.mesh.coord(:,1),T.mesh.coord(:,2),T.mesh.coord(:,3),...
            uZ, 'FaceColor','interp','EdgeColor','none');
colorbar
colormap turbo
xlabel('x')
ylabel('y')
zlabel('z')
axis('equal')
view(45,30)
title('displacement uZ(m)')
