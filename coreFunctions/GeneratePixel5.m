function [x,Density,sizeOfLattice] = GeneratePixel5(n,address,radius)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% n is the number of pixel along each axis
% address is the file location of wireframe
% density is the relative density
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin <4
    dof = 2;
end
if isscalar(n)
    n = n*ones(1,dof);
else
    dof = ndims(n);
end

%% Get truss topology
[node,connectivity] = ReadTrussTopology(address,dof);
mins = min(node,[],1);maxs = max(node,[],1);

%% According to size of lattice and number of elements get mesh
ds = (maxs-mins)./n;                               % initial size of voxels
[pixelNodes,pixelConn] = getMesh(dof,n,ds);

%% Allowing distinct radii for each truss member
if isscalar(radius), radius = radius*ones(size(connectivity,1),1);
elseif numel(radius) ~= size(connectivity,1), disp('Number of truss radii can either be a scalar or equal to number of truss members');
end

%% Get nodes inside topology + radius -> Get full pixels
[pixelNodeIn,signedDist] = getPtsInPeriodic(pixelNodes,radius,connectivity,node);
fullPixel = all(pixelNodeIn(pixelConn)== 1,2);
partialPixel = any(pixelNodeIn(pixelConn)== 0,2)&~fullPixel;

%% Find pts of 0 level contour and element boundaries
C = contourc(pixelNodes(1:(n(1)+1),1),pixelNodes(1:(n(1)+1):end,2),reshape(signedDist,n(1)+1,n(2)+1)',[0,0]);

%% Find partial density from contour data
if ~isempty(C)
    % Get intersection (bet contour and mesh connectivity)
    intersectionPts = zeros(size(C'));
    ctr = 1;ctr2 = 0;
    while(1)
        intersectionPts(ctr2+(1:C(2,ctr)),:) = C(:,ctr+(1:C(2,ctr)))';
        ctr2 = ctr2+C(2,ctr);
        ctr = ctr+C(2,ctr)+1;
        if ctr >size(C,2)
            break;
        end
    end
    intersectionPts(ctr2+1:end,:) = [];
    xs = intersectionPts(:,1);ys = intersectionPts(:,2);
    
    %Get intersection pt as multiple of element size
    fracElemX = xs/ds(1);fracElemY = ys/ds(2);
    
    %Fix rounding off error
    fracElemX(fracElemX<0) = 0;fracElemY(fracElemY<0) = 0;
    fracElemX(fracElemX>n(1)) = n(1);fracElemY(fracElemY>n(2)) = n(2);
    
    %Find intersection pts. lying on the elem boundaries
    tol = 1e-6;
    elemBoundX = checkEquality(round(fracElemX),(fracElemX),tol);
    elemBoundY = checkEquality(round(fracElemY),(fracElemY),tol);
    
    %Find elem indices
    is = [ceil(fracElemX(~elemBoundX&~elemBoundY));     %Those within elems
        fracElemX(elemBoundX);                          %Those on x nodes
        fracElemX(elemBoundX)+1;
        repmat(ceil(fracElemX(elemBoundY)),2,1)];       %Those on y nodes
    js = [ceil(fracElemY(~elemBoundX&~elemBoundY));     %Those within elems
        repmat(ceil(fracElemY(elemBoundX)),2,1);        %Those on x nodes
        fracElemY(elemBoundY);                          %Those on y nodes
        fracElemY(elemBoundY)+1];
    interNodes = [intersectionPts(~elemBoundX&~elemBoundY,:);
        repmat(intersectionPts(elemBoundX,:),2,1);
        repmat(intersectionPts(elemBoundY,:),2,1)];
    
    %Eliminate element indices not possible
    nonExistent = is<1|is>n(1)|js<1|js>n(2);
    is1 = round(is(~nonExistent));  js1 = round(js(~nonExistent));
    interNodes1 = interNodes(~nonExistent,:);
    
    %Get corresponding elements
    ijElems = is1+(js1-1)*n(1);
    
    %Remove element which are fully filled or fully empty
    partialElems = find(partialPixel(:));
    iA = ismember(ijElems,partialElems);
    ijElems1 = ijElems(iA);
    interNodes2 = interNodes1(iA,:);
    
    %Assign densities using convhull function to find area
    densityElem = zeros(size(fullPixel(:),1),1);
    densityElem(fullPixel) = 1;
    for el = 1:size(densityElem,1)
        if ~partialPixel(el), continue; end
        thisElemIntersNodes = interNodes2(ijElems1 == el,:);
        if isempty(thisElemIntersNodes), continue; end
        thisElemNodes = pixelConn(el,:);thisElemNodesIn = pixelNodeIn(pixelConn(el,:));
        thisElemNodeVal = pixelNodes(thisElemNodes,:);
        try  [~,ar] =convhull([thisElemNodeVal(thisElemNodesIn,:);thisElemIntersNodes]);
            densityElem(el) = ar/prod(ds);
        catch
            %When there are only two nodes
            %  disp('issue');
        end
    end
    x = reshape(densityElem,n);
else
    x = double(reshape(fullPixel,n));
end

Density = sum(x,'all')/prod(n);% calculate the relative density
sizeOfLattice = range(node,1);
x = permute(x,[2,1,3]);%To get more intuitive orientation (y-x-z) of pixels/voxels
x = x(end:-1:1,:,:);%Since the elements needs to stack from bottom to up and matrix is reshaped top to bottom

%------------------------------------------------------%
function [pixelNodes,pixelConn] = getMesh(dof,n,scale)
% Nodes are ordered as x -> y -> z
pixelNodes = zeros(prod(n+1),dof);
for i = 1:dof
    pixelNodes(:,i) = repmat(kron(0:n(i),ones(1,prod(n(1:i-1)+1))),1,prod(n(i+1:end)+1))';
end
pixelNodes = pixelNodes.*scale;
nodenrs = reshape(1:prod(n+1),n+1);
if dof == 2
    edofVec = reshape(nodenrs(1:end-1,1:end-1),prod(n),1);
    pixelConn = edofVec + [0 1 n(1)+2 n(1)+1];
elseif dof == 3
    edofVec = reshape(nodenrs(1:end-1,1:end-1,1:end-1),nel,1);
    pixelConn = edofVec + reshape([0 1 n(1)+2 n(1)+1]'+[0 prod(n(1:2)+1)],1,[]);
end
%------------------------------------------------------%
function ret = checkEquality(a,b,tol)
ret = and(a>b-tol,a<b+tol);
