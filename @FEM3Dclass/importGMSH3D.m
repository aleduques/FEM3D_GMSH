function [T, ComplementaryInformation] = importGMSH3D(filename)
%
% load the msh from GMSH,
% (format 4.1) specified in filename
% and set the mesh
%
% Outputs:
% ComplementaryInformation
% T = obj.msh
%
% Fields:
%
%  nNodes
%  deg
%  coord
%  trB
%  domBd
%  ttrh
%  domain
%  faces
%  ttrh2faces
%  faces2ttrh
%  trBNormal
%  facenormal
%  detBk
%


ComplementaryInformation = containers.Map('KeyType','char','ValueType','any');
meshFile =  fileread(filename);
meshFile =  regexp(meshFile, '\n', 'split');

% T : mesh struct
T.nNodes  = [];
T.deg     = [];
T.coord   = [];
T.trB     = [];
T.domBd   = [];
T.ttrh    = [];
T.domain  = [];
T.faces   = [];
T.ttrh2faces  = [];
T.faces2ttrh  = [];
T.nNodes   = [];
T.trBNormal   = [];
T.facenormal  = [];
T.detBk    = [];



%%%%%%%%%% We will check if everything is here:

disp('Checking Format')
if  isnan(strfind(meshFile{1},'$MeshFormat'))
    error('File format is not GMSH')
end
if  isnan(strfind(meshFile{2},'4.1'))
    warning('Mesh format is not apparently 4.1.')
end


% First thing: store the lines with the main information
sEntities      = findStringv2(meshFile,'$Entities',3);
eEntities      = findStringv2(meshFile,'$EndEntities',3);
sPhysicalNames = findStringv2(meshFile,'$PhysicalNames',3);
ePhysicalNames = findStringv2(meshFile,'$EndPhysicalNames',3);
sNodes         = findStringv2(meshFile,'$Nodes',3);
eNodes         = findStringv2(meshFile,'$EndNodes',3);
sElements      = findStringv2(meshFile,'$Elements',3);

if ~isnan(sPhysicalNames*ePhysicalNames)
    %% Insert Physical objects
    surfacesNum = []; surfacesPTags = {};
    volumesNum  = []; volumesPTags  = {};
    nPhysicalNames=str2num(meshFile{sPhysicalNames+1});
    for r = 1: nPhysicalNames
        aux    = regexp(meshFile{sPhysicalNames+1+r}, ' ', 'split');
        aux{1} = str2num(aux{1}); aux{2} = str2num(aux{2});
        tag = '';
        for i=3:length(aux)
            tag = [tag ' ' aux{i}];
        end
        if tag(end-1)=='"'
            tag = tag(3:end-2);
        else 
            tag = tag(3:end-1);
        end
        ComplementaryInformation(tag)=[];
        if aux{1} == 2
            surfacesNum(end+1)   = aux{2};
            surfacesPTags{end+1} = tag;
        elseif aux{1} == 3
            volumesNum(end+1)   = aux{2};
            volumesPTags{end+1} = tag;
        end
    end
end
%%%%% Entities
if ~isnan(sEntities*eEntities)
    aux = str2num(meshFile{sEntities+1});

    % Surfaces
    sSurfaces = sEntities+aux(1)+aux(2)+2;
    eSurfaces = sSurfaces + aux(3)-1;
    %Volumes
    sVolumes  = sEntities+aux(1)+aux(2)+aux(3)+2;
    eVolumes  = sVolumes+aux(4)-1;

    % Surfaces
    for i=sSurfaces:eSurfaces
        v = str2num(meshFile{i});
        nPhysicalTags = v(8);
        if  nPhysicalTags > 0
            for j =1:nPhysicalTags
                ComplementaryInformation(surfacesPTags{surfacesNum == v(8+j)}) = ...
                    [ComplementaryInformation(surfacesPTags{surfacesNum == v(8+j)}) v(1)];
            end
        end
    end
    % Volumes
    for i=sVolumes:eVolumes
        v = str2num(meshFile{i});
        nPhysicalTags = v(8);
        if  nPhysicalTags >0
            for j =1:nPhysicalTags
                ComplementaryInformation(volumesPTags{volumesNum == v(8+j)}) = ...
                    [ComplementaryInformation(volumesPTags{volumesNum == v(8+j)}) v(1)];
            end
        end
    end
end

if ~isnan(sNodes*eNodes)
    v = meshFile{sNodes+1};
    v = str2num(v);
    T.coord = nan(v(4),3);
    j = sNodes+2;
    while j<eNodes
        aux = str2num(meshFile{j});
        if aux(3) ~=0
            error('Parametric nodes are not supported')
        end
        % read the positions
        if aux(4)~=0
            indNodes = j+(1:aux(4));
            indCoord = indNodes(end)+(1:aux(4));

            indNodes = cellfun(@str2num, meshFile(indNodes));
            auxCoord = meshFile(indCoord); 
            str = strjoin(auxCoord,'\n');
            str = textscan(str,'%f');
            Coord = reshape(str{1}, [], numel(auxCoord))';
            T.coord(indNodes,:) = Coord;
        else
            warning(['Domain without proper nodes: ' num2str(aux)])
        end
        j = j+2*aux(4)+1;
    end
end

%%%%% Elements
disp('Checking Elements')
nBlocks = str2num(meshFile{sElements+1});
nBlocks = nBlocks(1);

pointerLine = sElements+2;
% supported Elements
supElements = [2  9  21   23 ;  % P1,P2, P3, P4 triangles
               4 11  29   30 ]; % P1, P2, P3, P4 tetrahedra
%Elements' fields initialization
T.trB   = []; T.domBd  = [];
T.ttrh  = []; T.domain = [];
T.faces = []; T.ttrh2faces = []; T.faces2ttrh = [];

for j=1:nBlocks
    aux = str2num(meshFile{pointerLine});
    deg = aux(3);
    supElementsTest = (supElements ==deg);
    indElements     = pointerLine+(1:aux(4));
    auxElm = meshFile(indElements); 
    str = strjoin(auxElm,'\n');
    str = textscan(str,'%f');
    elements = reshape(str{1}, [], numel(auxElm))';
    elements(:,1)   =[];
    if sum(supElementsTest(1,:))~=0  % 2D -> Face
        T.trB    = [T.trB; elements];
        T.domBd  = [T.domBd; ones(size(elements,1),1)*aux(2)];
    elseif sum(supElementsTest(2,:))~=0
        T.ttrh   = [T.ttrh; elements];
        T.domain = [T.domain; ones(size(elements,1),1)*aux(2)];
    end
    pointerLine  = pointerLine+aux(4)+1;
end


% Check if there are "ghost nodes"
nNodes     = max(max(T.ttrh));
nodesTtrh  = unique(T.ttrh(:));
nodesTrBd  = unique(T.trB(:));
if ~isempty(setdiff(nodesTrBd,nodesTtrh))
    error('Nodes in triangles which are not in ttrh')
end
[~, indghostNodes]= setdiff(1:nNodes,nodesTtrh);
if ~isempty(indghostNodes)
    % Fixing
    warning('ghost nodes... fixing')
    p      = nan(1,nNodes);
    p(nodesTtrh) = 1:length(nodesTtrh);
    T.coord(indghostNodes,:)=[];
    T.ttrh = p(T.ttrh);
    T.trB  = p(T.trB);
    nNodes = max(max(T.ttrh));
end
T.nNodes = nNodes;
switch(size(T.ttrh,2))
    case(4)
        T.deg =1;
    case(10)
        T.deg =2;
    case(20)
        T.deg =3;
    case(35)
        T.deg =4;
end



%%% Boundary faces
v21 = T.coord(T.trB(:,1),:) - T.coord(T.trB(:,2),:);
v31 = T.coord(T.trB(:,1),:) - T.coord(T.trB(:,3),:);
T.trBNormal = cross(v21, v31, 2);

%%% Tetrahedron faces

nElement = size(T.ttrh,1);
indAux   = [2 3 4; 1 4 3; 1 2 4; 1 3 2]; 
T.faces  = [T.ttrh(:,indAux(1,:)); ...
            T.ttrh(:,indAux(2,:)); ...
            T.ttrh(:,indAux(3,:)); ...
            T.ttrh(:,indAux(4,:))];
% We use only the faces' vertices
facesAux = sort(T.faces(:,1:3),2);     
[~,p1, p3]  = unique(facesAux,'rows', 'first');
[~,p2]      = unique(facesAux,'rows', 'last');
T.faces     = T.faces(p1,:);  
ind = zeros(nElement,4); ind(:) = 1:numel(T.ttrh(:,1:4));
T.ttrh2faces  = p3(ind);
p1  = mod(p1,nElement)    ; p1(p1==0) = nElement;   
p2  = mod(p2,nElement)    ; p2(p2==0) = nElement;  
T.faces2ttrh = [p1, p2];
ind = p1==p2;
T.faces2ttrh(ind,2) = nan; 


v21 = T.coord(T.faces(:,1),:) - T.coord(T.faces(:,2),:);
v31 = T.coord(T.faces(:,1),:) - T.coord(T.faces(:,3),:);
T.facenormal = cross(v21, v31, 2);

%%% Elements

% edges' vectors
v12 = T.coord(T.ttrh(:,2),:)-T.coord(T.ttrh(:,1),:);
v13 = T.coord(T.ttrh(:,3),:)-T.coord(T.ttrh(:,1),:);
v14 = T.coord(T.ttrh(:,4),:)-T.coord(T.ttrh(:,1),:);

% Determinant
T.detBk = dot(v12, cross(v13, v14, 2),2);
end

function [goAhead] = findStringv2(f,str,stLine)
    goAhead = stLine;
    while isempty(strfind(f{goAhead},str)) && goAhead<length(f)
        goAhead = goAhead+1;
    end
    if goAhead ==length(f) && isempty(strfind(f{goAhead},str))
        goAhead = nan;
    end
end
