T = FEM3Dclass('CylinderFins.msh');
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
        T.mesh.coord(:,3),'EdgeColor','k')
trimesh(ConvectSur, T.mesh.coord(:,1), T.mesh.coord(:,2), ...
        T.mesh.coord(:,3),'EdgeColor','k')
view(3)
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

disp('System assembly')
M = femMassMatrix(T, @(x,y,z) 1+x*0); M = rho*Cp*M;
S = femStressMatrix(T); S = k*S;
[t,MR] = femRobin(T, robin, gR, 'alpha',alpha); t = t(inD);
A  = M+0.5*tau*(S+MR); AnD =  A(inD,inD);
A2 = M-0.5*tau*(S+MR); 
disp('done')

% Can be used to decrease time computation of the system
% [l,d,p] = ldl(AnD);
clear  M S

maxDiff = [];
checktimes = [1 8 20 40];
for tn1 = t0+tau:tau:tf
    r   = A2(inD,:)*Uend+0.5*tau*(t+t)-A(inD,iD)*Uend(iD);
    Uend2 = AnD\r;
    maxDiff = [maxDiff max(abs(Uend(inD)-Uend2))];
    Uend(inD) = Uend2;
    if any(abs(tn1- checktimes)<1e-14)
        figure;
        trisurf(T.mesh.trB, T.mesh.coord(:,1)*Length,...
              T.mesh.coord(:,2)*Length,T.mesh.coord(:,3)*Length,...
              Uend-273, 'FaceColor','interp','EdgeColor','none');
        colormap turbo
        caxis([20 300])
        colorbar
        axis equal
        view(3)

        fig = figure;
        trisurf(T.mesh.trB, T.mesh.coord(:,1)*Length,...
              T.mesh.coord(:,2)*Length,T.mesh.coord(:,3)*Length,...
              Uend-273, 'FaceColor','interp','EdgeColor','none');
        colormap turbo
        caxis([20 300])
        colorbar
        axis equal
        view(0, 90)

        disp('press enter to continue')
        pause
    end
end

%%
figure
loglog(t0+tau:tau:tf, maxDiff, 'r-', 'linewidth',0.6)
ylabel('$\Delta T(^\circ C)$', 'interpreter','Latex')
xlabel('t')