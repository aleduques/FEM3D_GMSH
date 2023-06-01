function  b = femSourceTerm(obj, f)
% Returns the source vector of the PDE with source f
%
% Inputs:
% 
% -     f: Source function handle. f=f(x,y,z) or f=f(x,y,z,domain)
%
% Outputs:
%
% -     b: nNodesx1 source term.

% Geometry data
T = obj.mesh;
nTtrh  = size(T.ttrh,1);
nNodes = T.nNodes;
dofK = size(T.ttrh,2); 

% Cubature rule
degQuad      = max(2*T.deg-1,3);
quadRuleTtrh = obj.quadRule3D;
nodesQuad    = quadRuleTtrh{degQuad}.nodes; 
weights      = quadRuleTtrh{degQuad}.weights;  
nQ           = length(weights);

% Evaluate Interpolation functions
PkValues = zeros(dofK, nQ);
for i = 1:length(obj.Nj3D)
    PkValues(i,:) = obj.Nj3D{i}(nodesQuad(1,:),nodesQuad(2,:), ...
                                nodesQuad(3,:),nodesQuad(4,:)); 
end

% Cubature nodes in physical domain
px = T.coord(:,1); py = T.coord(:,2); pz = T.coord(:,3);
nodesX = px(T.ttrh(:,1:4))*nodesQuad;
nodesY = py(T.ttrh(:,1:4))*nodesQuad;
nodesZ = pz(T.ttrh(:,1:4))*nodesQuad;


if nargin(f) == 3
    val = f(nodesX, nodesY, nodesZ); 
elseif nargin(f) == 4
    val = f(nodesX, nodesY, nodesZ, T.domain*ones(1,nQ));
end
val = val.*T.detBk;

% Cubature Formula
PkValues = repmat(PkValues,nTtrh,1); 
val = repelem(val,dofK,1); 
val = val.*PkValues;
val = val*weights(:)/6; 

% Source vector assembly
indi = T.ttrh'; 
indi = indi(:); 
b =  accumarray(indi, val, [nNodes, 1]);
end

