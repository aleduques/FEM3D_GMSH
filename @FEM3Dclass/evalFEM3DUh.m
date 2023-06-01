function [val,Meval,barPt,indTtrh,indPtError] = evalFEM3DUh(uh,obj,pt,varargin)
%
% evaluate the FE solution uh defined on the mesh obj.mesh at pt
%
% [val,Meval,barPt,indTtrh,indPtError] = evalFEM3DUh(uh,obj,pt,list)
%
% input:
%        uh:  FE
%       obj:  FEM3Dclass
%        pt:  n x 3 points or nx4 where uh is going to be evaluated   
%      list:  tetrahedra where the points belongs to
%             (optional; for speeding up the computations)
%
% output:
%       val:  value of uh
%     Meval:  matrix to perform such a evaluation (via val = M*uh)
%     barpt:  The barycentric coordinates 4 x nPt
%   indTtrh:  Element(s) the point(s) belong to
% indPtError: Points for which not tetrahedra containing it was found


if isempty(uh)
    uh = obj.mesh.nNodes*nan;
end

if nargin <4 || isempty(varargin{1})
    list = 1:size(obj.mesh.ttrh,1);
else
    list = varargin{1};
end

if isempty(obj.Nj3D)
    % Local basis must be loaded
    obj.basisNj;
end
list = unique(list);

indPtError = [];

%% Find Elements containing the points
x = obj.mesh.coord(:,1);
y = obj.mesh.coord(:,2);
z = obj.mesh.coord(:,3);

x = x(obj.mesh.ttrh(list,1:4).');
y = y(obj.mesh.ttrh(list,1:4).');
z = z(obj.mesh.ttrh(list,1:4).');

% Block diagonal matrix to compute the barycentric coordinates
nT = length(list);
jaux = 1:nT*4;
jaux = [1 1 1 1]'*jaux(:)';

iaux = 1:4:nT*4;
iaux = kron(iaux, [1 1 1 1]);
iaux= bsxfun(@plus,[0 1 2 3]',iaux); 


matrix = sparse( iaux(:), jaux(:),[ones(1,4*nT); x(:)'; y(:)';z(:)'], 4*nT,4*nT);



nPt   = size(pt,1);
bt = kron(ones(nT,1),[pt(:,1).^0  pt ]');

aux = matrix\bt;
ind = abs(aux(1:4:end,:))+abs(aux(2:4:end,:)) + ...
    abs(aux(3:4:end,:))+abs(aux(4:4:end,:));

[aux2,indTtrh] = min(abs(ind-1));

Meval = sparse(nPt,obj.mesh.nNodes);
ind3= sub2ind(size(ind),indTtrh,1:nPt);

% barycentric coordinates
barPt= [aux(4*(ind3-1)+1); ...
        aux(4*(ind3-1)+2); ...
        aux(4*(ind3-1)+3); ...
        aux(4*(ind3-1)+4)];
% in tetrahedra indTtrh

indPtError = abs(aux2)>1e-8; 
if sum(indPtError)>0 & nargout<5
    warning('We haven''t found a tetrahedra for the following points')
    disp(pt(indPtError,:))
    disp('Using the closest tetrahedra instead')
    disp('with the following barycentric cooordinates')
    disp( barPt(:,indPtError)')
end

aux = zeros(size(obj.mesh.ttrh,2),nPt);
for j = 1:length(obj.Nj3D)
    aux(j,:)= obj.Nj3D{j}(barPt(1,:),barPt(2,:),...
                          barPt(3,:),barPt(4,:));
end
% index
indi = kron(1:nPt,ones(1,size(obj.mesh.ttrh,2))) ;
indj = obj.mesh.ttrh(list(indTtrh),:)'; indj = indj(:);
Meval = sparse(indi,indj,aux(:),nPt,obj.mesh.nNodes);

% Evaluation:

val = Meval*uh;

end

