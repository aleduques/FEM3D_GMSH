% T = FEM3Dclass('CylinderFins.msh');
T = FEM3Dclass('CylinderNoFins.msh');
%% Properties

% density 
rho = 7200;
% Thermal Conductivity
k = 52;
% Heat capacity
Cp = 447;


% Interior surface Temperature
TInt = 300; TInt = TInt+273;
% Convection
h = 5;
TInf = 27; TInf = TInf+273;
gR = @(x,y,z) h*TInf+x*0;
alpha = @(x,y,z) h+x*0;


%% Show geometry
% T.ComplementaryInformation.keys
FixedTempSur = T.ComplementaryInformation('Interior');
FixedTempSur = T.mesh.trB(ismember(T.mesh.domBd, FixedTempSur),:);

ConvectSur =  T.ComplementaryInformation('Exterior');
robin = ismember(T.mesh.domBd, ConvectSur);
ConvectSur = T.mesh.trB(ismember(T.mesh.domBd, ConvectSur),:);


figure
hold on
axis 'equal'
trimesh(FixedTempSur, T.mesh.coord(:,1), T.mesh.coord(:,2), ...
        T.mesh.coord(:,3),'FaceColor','red','EdgeColor','none')
trimesh(ConvectSur, T.mesh.coord(:,1), T.mesh.coord(:,2), ...
        T.mesh.coord(:,3),'FaceColor','blue','EdgeColor','none')
view(3)
%% Geometry Scaling
Length = 1000; %mm->m
T.mesh.coord = T.mesh.coord/Length;
T.mesh.detBk  = T.mesh.detBk/Length^3;
T.mesh.trBNormal  = T.mesh.trBNormal/Length^2;
%% Problem solving

S = femStressMatrix(T);
[t, MR] = femRobin(T, robin, gR, 'alpha',alpha);

A = k*S+MR;
b = t;
clear s MR t

u = zeros(T.mesh.nNodes,1);

% Dirichlet nodes
iD = unique(FixedTempSur(:));
inD = (1:T.mesh.nNodes)'; inD = setdiff(inD,iD);
u(iD) = TInt;


b = b(inD); b = b-A(inD,iD)*u(iD);
A = A(inD,inD);

u(inD) = A\b;

%% Solution representation
figure
trisurf(T.mesh.trB, T.mesh.coord(:,1)*Length,T.mesh.coord(:,2)*Length,T.mesh.coord(:,3)*Length,...
            u-273, 'FaceColor','interp','EdgeColor','none');
colormap turbo
colorbar
caxis([293,300])
axis equal
view(-23.5629,5.4)

figure
trisurf(T.mesh.trB, T.mesh.coord(:,1)*Length,T.mesh.coord(:,2)*Length,T.mesh.coord(:,3)*Length,...
            u-273, 'FaceColor','interp','EdgeColor','none');
colormap turbo
colorbar
caxis([293,300])
axis equal
view(0, 90)

