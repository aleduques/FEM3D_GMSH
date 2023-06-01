function [t,MR] = femRobin(obj, robin, gN, varargin)
% Returns the traction term and the boundary mass matrix of a 
% robin boundary condition of the form 
% -div(K grad u).n+ alpha*u=gN+alpha*u=gR.
%
% Inputs: 
% -      obj: object of the class FEM3Dclass.
% -      robin: logical vector indicating the elements on 
%        the robin boundary.
% -      gN: Scalar or vector function handle with the diffusion term of the
%        boundary condition. gN=gN(x,y,z) or gN=gN(x,y,z,domBd).
%           gN = @(...) gN or 
%           gN = @(...) [gNX gNY gNZ]
%
% Outputs:
% -      t: 1xnNodes traction term.
% -      MR: nNodesxnNodes boundary mass matrix.
%
%
% [t] = femRobin(obj, robin, gN). Constructs the traction term for the
% neumann problem -div(K*grad) u = gn. 
%
% [t, MR] =  femRobin(obj, robin, gR,'alpha', alpha). Construct the
%  traction term and the boundary mass matrix of the robin boundary
%  condition. alpha=alpha(x,y,z) or alpha=alpha(x,y,z,domBdR). 
%
%
% [t, MR] = femRobin(obj, robin, gN, 'alpha', alpha, 'exact', uExact).
% Construct the traction term with a known exact solution.  The exact
% solution to the BVP is uExact, and the Robin condition is given by
% gN+alpha*uExact=gR


p = inputParser;
addRequired(p,'obj')
addRequired(p,'robin')
addRequired(p, 'gN')
addParameter(p, 'alpha', @(x,y,z) 0.*x);
addParameter(p, 'exact', nan);
parse(p, obj, robin, gN, varargin{:});

% Obtention of robin parameters
alpha = p.Results.alpha;
u     = p.Results.exact;

% Geometry Data
T            = obj.mesh;
trBR         = T.trB(robin,:); 
domBdR       = T.domBd(robin);
ntrBR        = size(trBR,1);
normal       = T.trBNormal(robin,:);
detAk        = vecnorm(normal,2,2); 
nNodes       = T.nNodes;
dofA  = size(T.trB,2);


%%% Quadrature rule
degQuad    = max(2*T.deg-1,3);
quadRuleTr = obj.quadRule2D;
nodesQuad  = quadRuleTr{degQuad}.nodes;
weights    = quadRuleTr{degQuad}.weights; weights = weights(:);   
nQ         = length(weights);

% Evaluate Interpolation functions
PkValues = zeros(dofA, nQ);
for j = 1:length(obj.Nj2D)
    PkValues(j,:) = obj.Nj2D{j}(nodesQuad(1,:),...
                                nodesQuad(2,:),nodesQuad(3,:));  
end

% Quadrature nodes in physical domain
px = T.coord(:,1); py = T.coord(:,2); pz = T.coord(:,3);
nodesX = px(trBR(:,1:3))*nodesQuad;
nodesY = py(trBR(:,1:3))*nodesQuad;
nodesZ = pz(trBR(:,1:3))*nodesQuad;

%%%  Traction term
if nargin(gN) == 3
    val = gN(nodesX, nodesY, nodesZ); 
elseif nargin(gN) == 4
    val = gN(nodesX, nodesY, nodesZ, domBdR*ones(1,nQ)); 
end

if size(val,2)==nQ*3
    normal  = kron(normal, ones(1, size(val,2)/3));
    val     = val.*normal;
    val     = val(:,1:nQ)+val(:,nQ+1:2*nQ)+val(:,2*nQ+1:end);
elseif size(val,2)==nQ
    val     = val.*detAk;
end
if isa(u,'function_handle')
    if nargin(p.Results.exact)==3
        val = val + alpha(nodesX, nodesY, nodesZ).*...
                        u(nodesX, nodesY, nodesZ).*detAk;
    elseif nargin(p.Results.exact)==4
        val = val + alpha(nodesX, nodesY, nodesZ, domBdR).*...
                        u(nodesX, nodesY, nodesZ, domBdR).*detAk;
    end 
end
val = repmat(PkValues,ntrBR,1).*repelem(val,dofA, 1);
val = val*weights/2;

% Traction vector assembly
id  = trBR';
t   = accumarray(id(:), val, [nNodes, 1]);
clear val id normal

%%% Mass Boundary matrix
% Product of interpolation function at every quadrature node
PkProd = repmat(PkValues, dofA, 1).*repelem(PkValues, dofA, 1);
PkProd = repmat(PkProd, ntrBR, 1); 


if nargin(alpha) == 3
    val = alpha(nodesX, nodesY, nodesZ); 
elseif nargin(alpha) == 4
    val = alpha(nodesX, nodesY, nodesZ, domBdR*ones(1,nQ));
end
val = val.*detAk;

val = repelem(val, dofA^2, 1);
val = val.*PkProd; 
val = val*weights/2;

% Boundary mass matrix assembly
indi = repmat(trBR, 1, dofA);  indi=indi';     
indj = repelem(trBR, 1, dofA); indj=indj'; 
MR   = sparse(indi(:), indj(:), val(:), nNodes, nNodes);
end