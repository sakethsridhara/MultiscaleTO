%------------------ RegularMesher  ---------------------%
% Create arbitrary geometry using quad/triangular mesh
%-------------------------------------------------------------------------%
function [Node,Element,Supp,Load,P,ds] = RegularMesher(Domain,dim,NElem,AR)
BdBox = Domain('BdBox',dim);
w = BdBox(2)-BdBox(1);h = BdBox(4)-BdBox(3);
Area = w*h;
l = sqrt(Area/NElem/AR);
dx = l*AR;
dy = l;

nelx = 2*round(w/(2*dx));
nely = 2*round(h/(2*dy));
dx = w/nelx;
dy = h/nely;

nodex = (BdBox(1)+BdBox(2))/2+repelem((-dx*nelx/2:dx:dx*nelx/2)',nely+1,1);
nodey = (BdBox(3)+BdBox(4))/2+repmat((-dy*nely/2:dy:dy*nely/2)',nelx+1,1);
nodes = [nodex,nodey];
nodenrs = reshape(1:(1+nelx)*(1+nely),1+nely,1+nelx);

elemVec = reshape(nodenrs(1:end-1,1:end-1),nelx*nely,1);
elemMat = repmat(elemVec,1,4)+repmat([0 nely+[1 2] 1],nelx*nely,1);
elemCenter = [mean(nodex(elemMat),2),mean(nodey(elemMat),2)];
d = Domain('Dist',dim,elemCenter);
Element = elemMat(d(:,end)<0,:);
allNodesUnique = unique(Element(:));
Node = nodes(allNodesUnique,:);
elemMap(allNodesUnique) = 1:length(allNodesUnique);
Element = elemMap(Element);
P = elemCenter(d(:,end)<0,:);
ds = [dx,dy];
BC=Domain('BC',dim,{Node,Element}); Supp=BC{1}; Load=BC{2}; %Recover BC arrays
