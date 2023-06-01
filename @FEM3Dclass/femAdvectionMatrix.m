function [A] = femAdvectionMatrix(obj, v)
%  Returns the advection matrix of the PDE with fluid flow v
%
% Inputs:
% 
% -     v: Cell of function handles. 
%          v = {vx(x,y,z), vy(x,y,z), vz(x,y,z)} or
%          v = {vx(x,y,z,domain), vy(x,y,z,domain), vz(x,y,z,domain)}
%
% Outputs:
%
% -     A: nNodesxnNodes Advection Matrix.

% Geometry data
T = obj.mesh;
dofK = size(T.ttrh,2);
nNodes = T.nNodes;
nttrh=length(T.ttrh);

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

% Evaluate Interpolation functions
PkValues = zeros(dofK, nQ);
for j = 1:length(obj.Nj3D)
    PkValues(j,:) = ...
            obj.Nj3D{j}(nodesQuad(1,:),nodesQuad(2,:), ...
            nodesQuad(3,:),nodesQuad(4,:));
end

% Bk^-1.(T.detBk)
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
BsInv = sparse(rowIndex, colIndex, Bkp);


% Cubature nodes in physical domain
px = T.coord(:,1); py = T.coord(:,2); pz = T.coord(:,3);
nodesX = px(T.ttrh(:,1:4))*nodesQuad;
nodesY = py(T.ttrh(:,1:4))*nodesQuad;
nodesZ = pz(T.ttrh(:,1:4))*nodesQuad;

checkNargin = cellfun(@nargin,v);
if sum(sum(checkNargin~=checkNargin(1)))
    error(['The components of the velocity vector v should have' ...
        ' the same number of inputs'])
elseif nargin(v{1})==3
    vValx = v{1}(nodesX, nodesY, nodesZ);
    vValy = v{2}(nodesX, nodesY, nodesZ);
    vValz = v{3}(nodesX, nodesY, nodesZ);
elseif nargin(v{1})==4
    vValx = v{1}(nodesX, nodesY, nodesZ, T.domain);
    vValy = v{2}(nodesX, nodesY, nodesZ, T.domain);
    vValz = v{3}(nodesX, nodesY, nodesZ, T.domain);
end
vVal = [vValx vValy vValz];

%Cubature Formula
Aval = 0;
for i=1:nQ
    vnQ = sparse( indices, [ind ind ind], vVal(:, i:nQ:end), ...
                                               3*nttrh, 3*nttrh); 
    vnQ = BsInv*vnQ;
    vnQ = vnQ(:,ind);
    
    a1 = full(diag(vnQ(indices(:,1),:)));
    a2 = full(diag(vnQ(indices(:,2),:)));
    a3 = full(diag(vnQ(indices(:,3),:)));

    Aval = Aval+ 1/6*weights(i)*( ...
                 kron(a1', PkValues(:,i)*Sx(:,i)') + ...
                 kron(a2', PkValues(:,i)*Sy(:,i)') + ...
                 kron(a3', PkValues(:,i)*Sz(:,i)'));
end
clear vnQ a1 a2 a3

% Advection matrix assembly
[j, i]  = meshgrid(1:dofK,1:dofK);
indi    = zeros(dofK, dofK*nttrh); indj = indi;
indi(:) = T.ttrh(:,i)'; 
indj(:) = T.ttrh(:,j)';
A       = sparse(indi, indj, Aval, nNodes,nNodes); 
end

