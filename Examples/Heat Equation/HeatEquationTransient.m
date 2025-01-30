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
tic
T = FEM3Dclass('cylinderFins.msh');
toc
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

% Initial condition
T0 = 20; T0 = T0+273;
% time step
tau = 0.05;
t0 = 0;
tf = 120;


%% Show geometry

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

%set(gca,'XColor', 'none','YColor','none', 'ZColor','none')
view(120,10)
exportgraphics(gcf,"finned_cylinder.png", 'Resolution', 600);
%% Geometry Scaling
Length = 1000; %mm->m
T.mesh.coord = T.mesh.coord/Length;
T.mesh.detBk  = T.mesh.detBk/Length^3;
T.mesh.trBNormal  = T.mesh.trBNormal/Length^2;
 
%%
iD = unique(FixedTempSur(:));
inD = (1:T.mesh.nNodes)'; inD = setdiff(inD,iD);
Uend = zeros(T.mesh.nNodes, 1);


% Initial conditions
Uend = Uend+T0;
% Dirichlet conditions
Uend(iD) = TInt; 

%% Compute matrices
disp('System assembly')
tic
M = femMassMatrix(T, @(x,y,z) 1+x*0); M = rho*Cp*M;
toc
S = femStressMatrix(T); S = k*S;
[t,MR] = femRobin(T, robin, gR, 'alpha',alpha); t = t(inD);
disp("done")
%% Solve steady state
uSteady =  zeros(T.mesh.nNodes, 1);
uSteady(iD)=TInt;

A = S+MR;
b = t;


disp('Solving Steady state System')
b = b-A(inD,iD)*uSteady(iD);
A = A(inD,inD);
uSteady(inD) = A\b;
disp('Done')


%%
A  = M+0.5*tau*(S+MR); AnD =  A(inD,inD);
A2 = M-0.5*tau*(S+MR);


maxDiff = [];
comparisonSteady = [];
checktimes = [40 70 100];
tF = 300;
tic
disp("Computing transient state")
for tn1 = t0+tau:tau:tF
    r = A2(inD,:)*Uend + 0.5*tau*(t+t) - A(inD,iD)*Uend(iD);
    Uend2 = AnD\r;
    maxDiff = [maxDiff max(abs(Uend(inD)-Uend2))];
    Uend(inD) = Uend2;

    % % Check if current time matches any in checktimes
    % if any(abs(tn1 - checktimes) < 1e-14)
    %     time_label = round(tn1); 
    %     figure
    %     colormap jet
    %     trisurf(T.mesh.trB, T.mesh.coord(:,1)*Length, T.mesh.coord(:,2)*Length, T.mesh.coord(:,3)*Length, ...
    %             Uend - 273, 'FaceColor', 'interp', 'EdgeColor', 'none');
    %     set(gca, 'XColor', 'none', 'YColor', 'none', 'ZColor', 'none')
    %     cb = colorbar;
    %     cb.FontName = "TimesNewRoman";
    %     cb.FontSize = 14;
    %     ylabel(cb, "T [ºC]")
    %     caxis([293,300])
    %     grid off
    %     axis equal
    %     view(0, 90)
    %     exportgraphics(gcf, sprintf("top_%d.png", time_label), 'Resolution', 600);
    % 
    %     figure
    %     colormap jet
    %     trisurf(T.mesh.trB, T.mesh.coord(:,1)*Length, T.mesh.coord(:,2)*Length, T.mesh.coord(:,3)*Length, ...
    %             Uend - 273, 'FaceColor', 'interp', 'EdgeColor', 'none');
    %     set(gca, 'XColor', 'none', 'YColor', 'none', 'ZColor', 'none')
    %     cb = colorbar;
    %     cb.FontName = "TimesNewRoman";
    %     cb.FontSize = 14;
    %     ylabel(cb, "T [ºC]")
    %     caxis([293,300])
    %     grid off
    %     axis equal
    %     view(120, 10)
    %     exportgraphics(gcf, sprintf("view_%d.png", time_label), 'Resolution', 600);
    % end
    % comparisonSteady = [comparisonSteady max(abs(Uend - uSteady))];
end
toc
disp("Done")

