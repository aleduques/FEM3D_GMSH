% Make sure Elasticity is the current folder.
p = pwd ;
cd ..
cd ..
% Add class folder to path
addpath(pwd)
% Return to Elasticity Folder
cd(p)
%%

tic
T = FEM3Dclass('beam.msh');
toc
%%
figure
trBVertices = T.mesh.trB(:,1:3);
trimesh(trBVertices, T.mesh.coord(:,1),T.mesh.coord(:,2), T.mesh.coord(:,3),...
         'FaceColor',[.7 .7 .7],"EdgeColor", "black")
axis equal
hold on
fixed = T.ComplementaryInformation('Dirichlet');
fixed = T.mesh.trB(ismember(T.mesh.domBd, fixed),:);
fixedVertices = fixed(:,1:3);
trimesh(fixedVertices, T.mesh.coord(:,1),T.mesh.coord(:,2), T.mesh.coord(:,3),...
        'FaceColor','red',"EdgeColor", "black")
set(gca,'XColor', 'none','YColor','none', 'ZColor','none')
view(40, -50)
exportgraphics(gcf, "beam.png", 'Resolution', 600);

%% Geometry Scaling
Length = 1000; %mm->m
T.mesh.coord = T.mesh.coord/Length;
T.mesh.detBk  = T.mesh.detBk/Length^3;
T.mesh.trBNormal  = T.mesh.trBNormal/Length^2;
 
%% Nodos Dirichlet

nNodes = T.mesh.nNodes;
fixedNodes = unique(fixed(:));
iD = [fixedNodes; fixedNodes+nNodes; fixedNodes+2*nNodes];
inD = (1:3*nNodes)'; inD = setdiff(inD, iD);

%% Propiedades de materiales 

E = 210e9; %N/m^2
nu = 0.3;
rho = 7850; %kg/m3

lambda = E*nu/(1+nu)/(1-2*nu);
mu = E/2/(1+nu);

%% Ensamblaje de sistema
L2Mu = @(x,y,z) lambda + 2*mu + 0*x;
Mu  = @(x,y,z) mu+0*x;
L   = @(x,y,z) lambda+0*x;
zero = @(x,y,z) 0*x;

disp('Assemble stiffness matrix')
tic
A11 = {
    L2Mu, zero, zero;
    zero, Mu,   zero;
    zero, zero, Mu
};
S11 = femStiffnessMatrix(T,A11);

A22 = {
    Mu, zero, zero;
    zero, L2Mu,   zero;
    zero, zero, Mu;
};
S22 = femStiffnessMatrix(T,A22);

A33 = {
    Mu, zero, zero;
    zero, Mu,   zero;
    zero, zero, L2Mu
};
S33 = femStiffnessMatrix(T,A33);

A12 = {
    zero, L, zero;
    Mu,   zero,zero;
    zero,zero,zero
};
S12 = femStiffnessMatrix(T,A12);

A13 = {
    zero, zero, L;
    zero, zero,   zero;
    Mu, zero, zero
};
S13 = femStiffnessMatrix(T,A13);

A23 = {
    zero, zero, zero;
    zero, zero, L;
    zero, Mu, zero
};
S23 = femStiffnessMatrix(T,A23);

matrix = [S11 S12 S13
          S12' S22 S23
          S13' S23' S33];
stiffness = matrix(inD, inD);
toc
disp("done")
disp("Assemble Mass matrix")
tic
M = femMassMatrix(T,@(x,y,c) rho+0*x);

mass = blkdiag(M,M,M);
mass = mass(inD, inD);
disp("Done")
toc
%%

nEigs = 5;
disp("Compute resonance frequencies")
tic
[ v, d] = eigs(stiffness, mass, nEigs,'smallestabs');
toc
d = diag(d);
omega = sqrt(d);
frequency = sqrt(d)/2/pi;
%%
trB = T.mesh.trB;
trB = [trB(:, [1 4 6]);
       trB(:, [4 2 5]);
       trB(:, [6 5 3]);
       trB(:, [4 5 6])];
filename = {"first_mode.png", "second_mode.png", "third_mode.png", "fourth_mode.png", "fifth_mode.png"};
for i=1:nEigs
    uMode1 = zeros(nNodes*3,1);
    uMode1(inD) = v(:,i); 
    uMode1X = uMode1(1:nNodes);
    uMode1Y = uMode1(nNodes+1:2*nNodes);
    uMode1Z = uMode1(2*nNodes+1:3*nNodes);
    
    uModule = vecnorm([uMode1X, uMode1Y, uMode1Z],2,2);
    figure
    fig = gcf;
    trisurf(trB, T.mesh.coord(:,1)+uMode1X, T.mesh.coord(:,2)+uMode1Y,...
        T.mesh.coord(:,3)+uMode1Z, uModule, 'FaceColor', 'interp','Edgecolor','none')
    axis equal
    ax = gca; 
    axPos = ax.Position;
    ax.Position = [axPos(1) - 0.05, axPos(2), axPos(3), axPos(4)];
    set(gca,'XColor', 'none','YColor','none', 'ZColor','none')
    cb = colorbar;
    cbPos = cb.Position;
    cb.Position = [cbPos(1) - 0.1, cbPos(2), cbPos(3), cbPos(4)];
    cb.FontName = "TimesNewRoman";
    cb.FontSize = 14;
    grid off
    ylabel(cb, "$\|u\|$ [m]","Interpreter","latex")
    
    exportgraphics(gcf, filename{i}, 'Resolution', 600);

end
