function  Mc = femMassMatrix(obj, c)
% Returns the Mass matrix of the PDE with linear reaction c*u
%
% Inputs:
% 
% -     c: Reaction function handle. c=c(x,y,z) or c=c(x,y,z,domain)
%
% Outputs:
%
% -     Mc: nNodesxnNodes Mass Matrix.

% Geometry Data
T           = obj.mesh; 
nTtrh       = size(T.ttrh,1);
nNodes      = T.nNodes;
dofK = size(T.ttrh,2); 

%%% Cubature rule
degQuad      = max(2*T.deg-1,3);
quadRuleTtrh = obj.quadRule3D;
nodesQuad    = quadRuleTtrh{degQuad}.nodes; 
weights      = quadRuleTtrh{degQuad}.weights;   
nQ           = length(weights);

% Evaluate Interpolation functions
PkValues =zeros(dofK, nQ);
for i = 1:length(obj.Nj3D)
    PkValues(i,:) = obj.Nj3D{i}(nodesQuad(1,:),nodesQuad(2,:), ...
                                nodesQuad(3,:),nodesQuad(4,:)); 
end
PkProd = repelem(PkValues,dofK,1).*repmat(PkValues,dofK,1);

% Cubature nodes in physical domain
px = T.coord(:,1); py = T.coord(:,2); pz = T.coord(:,3);
nodesX = px(T.ttrh(:,1:4))*nodesQuad;
nodesY = py(T.ttrh(:,1:4))*nodesQuad;
nodesZ = pz(T.ttrh(:,1:4))*nodesQuad;

if nargin(c) == 3
    val = c(nodesX, nodesY, nodesZ);
elseif nargin(c) == 4
    val = c(nodesX, nodesY, nodesZ, T.domain*ones(1,nQ));
end
val = val.*T.detBk;

% Cubature formula
PkProd = repmat(PkProd, nTtrh, 1); 
val    = repelem(val, dofK^2, 1); 
val    = val.*PkProd; 
val    = val*weights(:)/6; 

% Mass matrix assembly
indi = repmat(T.ttrh,1,dofK);  indi = indi';     
indj = repelem(T.ttrh,1,dofK); indj = indj'; 
Mc   = sparse(indi(:), indj(:), val(:), nNodes, nNodes);
end



