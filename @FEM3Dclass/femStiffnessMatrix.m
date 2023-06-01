function [S] = femStiffnessMatrix(obj,K)
% Returns the Stiffness matrix of the PDE with diffusion term
% -div(K grad u)
%
% Inputs:
% 
% -     K: Diffusion cell function handle. 
%    K = {K11(x,y,z), K12(x,y,z), K13(x,y,z)
%         K21(x,y,z), K22(x,y,z), K23(x,y,z)
%         K31(x,y,z), K32(x,y,z), K33(x,y,z)} or
%    K = {K11(x,y,z,domain), K12(x,y,z,domain), K13(x,y,z,domain)
%         K21(x,y,z,domain), K22(x,y,z,domain), K23(x,y,z,domain)
%         K31(x,y,z,domain), K32(x,y,z,domain), K33(x,y,z,domain)}
%
% Outputs:
%
% -     S: nNodesxnNodes Stiffness Matrix.


% Geometry data
T           = obj.mesh; 
nttrh       = size(T.ttrh,1);
nNodes      = T.nNodes;
dofK = size(T.ttrh,2); 

%Cubature Rule
degQuad      = max(2*T.deg-1,3);
quadRuleTtrh = obj.quadRule3D;
nodesQuad    = quadRuleTtrh{degQuad}.nodes; 
weights      = quadRuleTtrh{degQuad}.weights;  
nQ           = length(weights);

% Evaluate gradient of Interpolation function 
gradPhi  = zeros(3*dofK,nQ);
for j=1:3
    for i=1:dofK
        gradPhi((j-1)*dofK+(i-1)+1,:) = ...
        obj.gradNj3D{j,i}(nodesQuad(1,:), nodesQuad(2,:), ...
                          nodesQuad(3,:), nodesQuad(4,:));
    end
end
Sx = gradPhi(1:dofK,:);
Sy = gradPhi(dofK+1:2*dofK,:);
Sz = gradPhi(2*dofK+1:end,:);
clear gradPhi

% Bk^-1.sqrt(T.detBk)
v12 = T.coord(T.ttrh(:,2),:)-T.coord(T.ttrh(:,1),:);
v13 = T.coord(T.ttrh(:,3),:)-T.coord(T.ttrh(:,1),:);
v14 = T.coord(T.ttrh(:,4),:)-T.coord(T.ttrh(:,1),:);

b1 = cross(v13,v14,2);
b2 = cross(v14,v12,2);
b3 = cross(v12,v13,2);

ind = 1:nttrh;
ind32 = 3*ind(:)-2;
ind31 = 3*ind(:)-1;
ind30 = 3*ind(:);
indices = [ind32 ind31 ind30];
rowIndex = repmat(indices(:),3,1);
colIndex = reshape(repmat(indices,3,1),[],1);

Bkp   = [b1; b2; b3];
Bkp   = Bkp./sqrt(repmat(T.detBk,3,1));
BsInv = sparse(rowIndex, colIndex, Bkp);

% Cubature nodes in physical domain
px = T.coord(:,1); py = T.coord(:,2); pz = T.coord(:,3);
nodesX = px(T.ttrh(:,1:4))*nodesQuad;
nodesY = py(T.ttrh(:,1:4))*nodesQuad;
nodesZ = pz(T.ttrh(:,1:4))*nodesQuad;

checkNargin = cellfun(@nargin,K);
if sum(sum(checkNargin~=checkNargin(1)))
    error(['The components of the diffusion matrix K should have' ...
        ' the same number of inputs'])
elseif nargin(K{1,1})==3
    K11 = K{1,1}(nodesX, nodesY, nodesZ);
    K12 = K{1,2}(nodesX, nodesY, nodesZ);
    K13 = K{1,3}(nodesX, nodesY, nodesZ);
    K21 = K{2,1}(nodesX, nodesY, nodesZ);
    K22 = K{2,2}(nodesX, nodesY, nodesZ);
    K23 = K{2,3}(nodesX, nodesY, nodesZ);
    K31 = K{3,1}(nodesX, nodesY, nodesZ);
    K32 = K{3,2}(nodesX, nodesY, nodesZ);
    K33 = K{3,3}(nodesX, nodesY, nodesZ);
elseif nargin(K{1,1})==4
    K11 = K{1,1}(nodesX, nodesY, nodesZ,T.domain);
    K12 = K{1,2}(nodesX, nodesY, nodesZ,T.domain);
    K13 = K{1,3}(nodesX, nodesY, nodesZ,T.domain);
    K21 = K{2,1}(nodesX, nodesY, nodesZ,T.domain);
    K22 = K{2,2}(nodesX, nodesY, nodesZ,T.domain);
    K23 = K{2,3}(nodesX, nodesY, nodesZ,T.domain);
    K31 = K{3,1}(nodesX, nodesY, nodesZ,T.domain);
    K32 = K{3,2}(nodesX, nodesY, nodesZ,T.domain);
    K33 = K{3,3}(nodesX, nodesY, nodesZ,T.domain);
end
clear px py pz nodesX nodesY nodesZ Bkp

% Cubature formula
Sval = 0;
for i=1:nQ
    KK   = [K11(:,i) K12(:,i) K13(:,i); ...
            K21(:,i) K22(:,i) K23(:,i); ...
            K31(:,i) K32(:,i) K33(:,i)];
      
    Cs   = sparse(rowIndex, colIndex, KK);
    Cs   = BsInv*Cs*BsInv';
    aux  = diag(Cs); 
    C11  = full(aux(1:3:end));
    C22  = full(aux(2:3:end));
    C33  = full(aux(3:3:end));
    aux  = diag(Cs,1);
    C12  = full(aux(1:3:end));
    C23  = full(aux(2:3:end));
    aux  = diag(Cs,-1);
    C21  = full(aux(1:3:end));
    C32  = full(aux(2:3:end));
    aux  = diag(Cs,2);
    C13  = full(aux(1:3:end));
    aux  = diag(Cs,-2);
    C31  = full(aux(1:3:end));
    
    Sval = Sval +  1/6*weights(i)*(   ...
        kron(C11',Sx(:,i)*Sx(:,i)') + ...
        kron(C22',Sy(:,i)*Sy(:,i)') + ...
        kron(C33',Sz(:,i)*Sz(:,i)') + ...
        kron(C12',Sx(:,i)*Sy(:,i)') + ...
        kron(C13',Sx(:,i)*Sz(:,i)') + ...
        kron(C21',Sy(:,i)*Sx(:,i)') + ...
        kron(C23',Sy(:,i)*Sz(:,i)') + ...
        kron(C31',Sz(:,i)*Sx(:,i)') + ...
        kron(C32',Sz(:,i)*Sy(:,i)')   ...
    );
end
clear C11 C22 C33 C12 C23 C21 C32 C13 C31 aux KK K11 K12 K13...
      K21 K22 K23 K31 K32 K33

% Stiffness matrix Assembly
[j, i]  = meshgrid(1:dofK,1:dofK);
indi    = zeros(dofK, dofK*nttrh); indj = indi;
indi(:) = T.ttrh(:,i)'; 
indj(:) = T.ttrh(:,j)';
S       = sparse(indi, indj, Sval, nNodes, nNodes);
end