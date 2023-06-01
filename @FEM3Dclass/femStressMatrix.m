function  S = femStressMatrix(obj)
% Returns the Stiffness matrix of the PDE with diffusion term
% -lap u
%
% Outputs:
%
% -     S: nNodesxnNodes Stiffness Matrix.

% Geometry Data
T           = obj.mesh; 
nttrh       = size(T.ttrh,1);
nNodes      = T.nNodes;
dofK = size(T.ttrh,2); 

[~,Sxx,Sxy,Sxz,Syy,Syz,Szz] = obj.local3DMatrices(T.deg);


%Ck matrix elements obtention
v12 = T.coord(T.ttrh(:,2),:)-T.coord(T.ttrh(:,1),:);
v13 = T.coord(T.ttrh(:,3),:)-T.coord(T.ttrh(:,1),:);
v14 = T.coord(T.ttrh(:,4),:)-T.coord(T.ttrh(:,1),:); 

c11 =  v13(:,3).^2.*(v14(:,1).^2+v14(:,2).^2)-        ...
       2.*v13(:,1).*v13(:,3).*v14(:,1).*v14(:,3)-2.*  ...
       v13(:,2).*v14(:,2).*(v13(:,1).*v14(:,1)+       ...
       v13(:,3).*v14(:,3))+v13(:,2).^2.*(v14(:,1).^2+ ...
       v14(:,3).^2)+v13(:,1).^2.*(v14(:,2).^2+        ...
       v14(:,3).^2); 
c11 = c11./T.detBk;
c22 =  v12(:,3).^2.*(v14(:,1).^2+v14(:,2).^2)-        ...
       2.*v12(:,1).*v12(:,3).*v14(:,1).*v14(:,3)-2.*  ...
       v12(:,2).*v14(:,2).*(v12(:,1).*v14(:,1)+       ...
       v12(:,3).*v14(:,3))+v12(:,2).^2.*(v14(:,1).^2+ ...
       v14(:,3).^2)+ v12(:,1).^2.*(v14(:,2).^2+       ...
       v14(:,3).^2); 
c22 = c22./T.detBk;
c33 =  v12(:,3).^2.*(v13(:,1).^2+v13(:,2).^2)-        ...
       2.*v12(:,1).*v12(:,3).*v13(:,1).*v13(:,3)-2.*  ...
       v12(:,2).*v13(:,2).*(v12(:,1).*v13(:,1)+       ...
       v12(:,3).*v13(:,3))+v12(:,2).^2.*(v13(:,1).^2+ ...
       v13(:,3).^2)+ v12(:,1).^2.*(v13(:,2).^2+       ...
       v13(:,3).^2); 
c33 = c33./T.detBk;
c12 =  -v12(:,3).*v13(:,3).*(v14(:,1).^2+v14(:,2).^2)+...
       v12(:,3).*(v13(:,1).*v14(:,1)+                 ...
       v13(:,2).*v14(:,2)).* v14(:,3)+                ...
       v12(:,2).*v14(:,2).*(v13(:,1).*v14(:,1)+       ...
       v13(:,3).*v14(:,3))-v12(:,2).*v13(:,2).*(      ...
       v14(:,1).^2+v14(:,3).^2)+                      ...
       v12(:,1).*(v13(:,2).*v14(:,1).*v14(:,2)+       ...
       v13(:,3).*v14(:,1).*v14(:,3)-                  ...
       v13(:,1).*(v14(:,2).^2+v14(:,3).^2)); 
c12 = c12./T.detBk;
c13 = v12(:,3).*(v13(:,1).*v13(:,3).*v14(:,1)+        ...
      v13(:,2).*v13(:,3).*v14(:,2)-                   ...
      v13(:,1).^2.*v14(:,3)- v13(:,2).^2.*v14(:,3))+  ...
      v12(:,1).*(-v13(:,2).^2.*v14(:,1)-              ...
      v13(:,3).^2.*v14(:,1)+                          ...
      v13(:,1).*v13(:,2).*v14(:,2)+                   ...
      v13(:,1).*v13(:,3).*v14(:,3))+                  ...
      v12(:,2).*(v13(:,1).*v13(:,2).*v14(:,1)-        ...
      v13(:,1).^2.*v14(:,2)+                          ...
      v13(:,3).*(-v13(:,3).*v14(:,2)+v13(:,2).*v14(:,3))); 
c13 = c13./T.detBk;
c23 = -v12(:,3).^2.*(v13(:,1).*v14(:,1)+              ...
      v13(:,2).*v14(:,2))+                            ...
      v12(:,1).*v12(:,3).*(v13(:,3).*v14(:,1)+        ...
      v13(:,1).*v14(:,3))+                            ...
      v12(:,2).*(v12(:,1).*v13(:,2).*v14(:,1)+        ...
      v12(:,1).*v13(:,1).*v14(:,2)+                   ...
      v12(:,3).*v13(:,3).*v14(:,2)+                   ...
      v12(:,3).*v13(:,2).*v14(:,3))-                  ...
      v12(:,2).^2.*(v13(:,1).*v14(:,1)+               ...
      v13(:,3).*v14(:,3))-                            ...
      v12(:,1).^2.*(v13(:,2).*v14(:,2)+               ...
      v13(:,3).*v14(:,3)); 
c23 = c23./T.detBk;

Sk = kron(c11', Sxx)+kron(c22', Syy)+kron(c33', Szz)+ ...
     kron(c12', Sxy+Sxy')+ kron(c13',Sxz+Sxz')+       ...
     kron(c23', Syz+Syz');


% Stiffness matrix assembly
[j,i] = meshgrid(1:dofK, 1:dofK);
indi  = zeros(dofK, dofK*nttrh); indj = indi;
indi(:) = T.ttrh(:,i)'; 
indj(:) = T.ttrh(:,j)';
S = sparse(indi, indj, Sk, nNodes, nNodes);
end
