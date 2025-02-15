%% 
%%% If plots are not showing in octave:
% graphics_toolkit ("gnuplot")

% Make sure Heat Equation is the current folder.
p = pwd ;
cd ..
cd ..
% Add class folder to path
addpath(pwd)
% Return to Heat Equation Folder
cd(p)

%%
T = FEM3Dclass('cylinderFins.msh');
% T = FEM3Dclass('cylinderNoFins.msh');
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
        T.mesh.coord(:,3),'FaceColor','red','EdgeColor','black')
trimesh(ConvectSur, T.mesh.coord(:,1), T.mesh.coord(:,2), ...
        T.mesh.coord(:,3),'FaceColor','#ECECEC','EdgeColor','black')

set(gca,'XColor', 'none','YColor','none', 'ZColor','none')
view(120,10)
%% Geometry Scaling
Length = 1000; %mm->m
T.mesh.coord = T.mesh.coord/Length;
T.mesh.detBk  = T.mesh.detBk/Length^3;
T.mesh.trBNormal  = T.mesh.trBNormal/Length^2;
%% Problem solving

disp('System Assembly')
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

disp('Solving System')
b = b(inD); b = b-A(inD,iD)*u(iD);
A = A(inD,inD);

u(inD) = A\b;
disp('done')
%% Solution representation
close all
figure
colormap jet
trisurf(T.mesh.trB, T.mesh.coord(:,1)*Length,T.mesh.coord(:,2)*Length,T.mesh.coord(:,3)*Length,...
            u-273, 'FaceColor','interp','EdgeColor','none');
set(gca,'XColor', 'none','YColor','none', 'ZColor','none')
cb = colorbar;
cb.FontName = "TimesNewRoman";
cb.FontSize = 14;
ylabel(cb, "T [ºC]")
caxis([293,300])
axis equal
grid off
view(0,90)
exportgraphics(gcf, "steady_top.png", 'Resolution', 600);

figure
colormap jet
trisurf(T.mesh.trB, T.mesh.coord(:,1)*Length,T.mesh.coord(:,2)*Length,T.mesh.coord(:,3)*Length,...
            u-273, 'FaceColor','interp','EdgeColor','none');
set(gca,'XColor', 'none','YColor','none', 'ZColor','none')
cb = colorbar;
cb.FontName = "TimesNewRoman";
cb.FontSize = 14;
ylabel(cb, "T [ºC]")
caxis([293,300])
axis equal
grid off
view(120,10)
exportgraphics(gcf, "steady_view.png", 'Resolution', 600);
