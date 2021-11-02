%% Types of MS
%1 = Plus, 2 = square, 3 = inverted Z; 4 = Z; 5 = N; 6 = inverted N; 7 =x;
% 8 = hexagon iso; 9 = rectangle orthotropic; 10 = reentrant hex; 11 = swastik;
% 12 = kagome iso; 13 = diamond orth; 14 = tri honeycomb almost iso;

%% Optimize
function OptimizeSD(dType, msType, deg, volFrac)
close all;
[dim,MSType,fileName,matProp,AR,nElems] = setProb(dType,msType);
%addpath('./PolyMesher/','../coreFunctions/');

[node,elem,supp,loads,~,ds] = RegularMesher(fileName,dim,nElems,AR);
[maxIter,objType,filterRad] = setOptim();
[densities] = runOptim(dim,MSType,deg,elem,node,ds,matProp,supp,loads,volFrac,maxIter,objType,filterRad);

%store the designVariables in a struct object
oldDesign.densities = densities;
oldDesign.node = node;
oldDesign.elem = elem;
oldDesign.ds = ds;
oldDesign.matProp = matProp;
oldDesign.MSType = MSType;
%Print design using microstructure topology and density info
printDesignNew1(dim,oldDesign);

%--------------------------------------------%
function [dim,MSType,file,matProp,AR,nElems]=setProb(dType,msType)
dim = 2;
E = 1; nu = 0.3; t = 1;
matProp = [1,E,nu,t];
MSType = msType;   % MS_type:
AR = getAspectRatio(MSType);
switch dType    %Set domain and nElems here
    case 1
        file = @CantileverDomain;   nElems = 5000;
    case 2
        file = @LBracketDomain;     nElems = 10000;
    case 3
        file = @MbbDomain;          nElems = 5000;
end
%--------------------------------------------%
function [maxIter,objType,filterRad] = setOptim()
%Set TO problem parameters
maxIter = 800;
objType = 1;                    %1 = compliance
filterRad = 1.2;                %filter radius wrt macro_elem size

%% Don't modify anything below this line
%--------------------------------------------%
function printDesignNew1(dim,oldDes)
folder = './latticeTopologies/';
filename = ['lattice2d_',num2str(oldDes.MSType)];
ext = '.txt';
address = [folder,filename,ext];
folder2 = './FitData/';
fileToWrite = 'thickVsDensity';
suffix = '.mat';
thicknessFile = [folder2 fileToWrite num2str(oldDes.MSType) suffix];
load(thicknessFile,'tfunc');
[nodeList] = ReadTrussTopology(address,dim);
mins = min(nodeList,[],1);maxs = max(nodeList,[],1);
scaling = (maxs - mins)./oldDes.ds;
radForAllElems = ppval(tfunc,oldDes.densities);
pElem = polyshape([mins(1) maxs(1) maxs(1) mins(1)],[mins(2) mins(2) maxs(2) maxs(2)]);
[polyshapeFunctionPos,polyshapeFunctionNeg] = createGradedMS_polyShapeFunc(oldDes.MSType);
psAll = repmat(polyshape(),size(oldDes.elem,1),numel(polyshapeFunctionPos));
psAllNeg = repmat(polyshape(),size(oldDes.elem,1),numel(polyshapeFunctionNeg));
psNegVoid = repmat(polyshape(),1,numel(polyshapeFunctionNeg));
psNegNonVoid = repmat(polyshape(),1,numel(polyshapeFunctionNeg));
for el = 1:size(oldDes.elem,1)
    if oldDes.densities(el)<0.01
        continue;
    elseif oldDes.densities(el)>0.99
        ps = pElem;
        psNeg = psNegVoid;
    else
        thisRad = mean(radForAllElems(el));
        ps = polyshapeFunctionPos(thisRad);
        if iscell(polyshapeFunctionNeg) 
            psNeg = psNegNonVoid;
            for i = 1:size(polyshapeFunctionNeg)
                psNeg(i) = polyshapeFunctionNeg{i}(thisRad);
            end
        else
            psNeg = polyshapeFunctionNeg(thisRad);
        end
    end
    psVec = scale(ps,1./scaling);
    transX = oldDes.node(oldDes.elem(el,1),1);
    transY = oldDes.node(oldDes.elem(el,1),2);
    psVec = translate(psVec,transX,transY);
    psAll(el,1:numel(psVec)) = psVec(:)';
    psVec = scale(psNeg,1./scaling);
    psVec = translate(psVec,transX,transY);
    psAllNeg(el,1:numel(psVec)) = psVec(:)';
end
%Get vertices of each polyshape object
var = arrayfun(@(x) {x.Vertices},psAll(:));
varNeg = arrayfun(@(x) {x.Vertices},psAllNeg(:));
%remove empty cells 
var = var(~cellfun('isempty', var));
varNeg = varNeg(~cellfun('isempty', varNeg));
%Create an array of vertices
vertices = cell2mat(var);
verticesNeg = cell2mat(varNeg);
%Maximum number of vertices in any polyshape object
maxPoly = max(arrayfun(@(x) size(x{1},1),var));
maxPolyNeg = max(arrayfun(@(x) size(x{1},1),varNeg));
%Vertex Serial number of polyshape object
nVerts = cumsum(arrayfun(@(x)size(x{1},1),var));
nVertsNeg = cumsum(arrayfun(@(x)size(x{1},1),varNeg));
%Create patches by writing vertex index of each patch no. start with 1
facesCell = cellfun(@(x) [1:size(x,1),1,NaN(1,maxPoly-size(x,1))],var,'UniformOutput',false);
facesCellNeg = cellfun(@(x) [1:size(x,1),1,NaN(1,maxPolyNeg-size(x,1))],varNeg,'UniformOutput',false);
%Update patch info by adding vertex serial number to vertex index
facesCell(2:end) = arrayfun(@(x,y) x{1}+y,facesCell(2:end),nVerts(1:end-1),'UniformOutput',false);
facesCellNeg(2:end) = arrayfun(@(x,y) x{1}+y,facesCellNeg(2:end),nVertsNeg(1:end-1),'UniformOutput',false);
%Get faces as matrix
faces = cell2mat(facesCell);
facesNeg = cell2mat(facesCellNeg);
figure('units','normalized','outerposition',[0 0 1 1]);
patch('faces',faces,'vertices',vertices,'faceColor','k','edgeColor','none','facealpha',1);
patch('faces',facesNeg,'vertices',verticesNeg,'faceColor','w','edgeColor','none','facealpha',1);
axis equal off tight;
drawnow;

function density = runOptim(dim,msType,deg,elem,node,ds,matProp,supp,loads,volFrac,maxIter,objType,filterRad)
nNode = size(node,1); nElem = size(elem,1);
sizeOfC = 3*(dim-1);    nIndDcomp = (sizeOfC*(sizeOfC+1))/2;
nDofs = dim*nNode;

dx = ds(1); dy = ds(2);

%Setting dim of each element and index vector for stiffness matrix K
elemDofs = reshape(reshape((2*reshape(elem',[],1)-[1,0])',[],1),8,[])';
iK = reshape(kron(elemDofs,ones(8,1))',64*nElem,1);
jK = reshape(kron(elemDofs,ones(1,8))',64*nElem,1);

%Computing element stiffness matrix template
map = {[1,1],[2,2],[3,3],[1,2],[1,3],[2,3]};
t = matProp(1,4);
[B,Det] = calc_Bdet1(elem,node,1);
Ke_template = zeros(8,8,nIndDcomp);
for i = 1:nIndDcomp
    C = zeros(sizeOfC);
    C(map{i}(1),map{i}(2)) = 1;C(map{i}(2),map{i}(1)) = 1;
    Ke_template(:,:,i) = Kelem_congruElem(C,t,B,Det{1,1});
end

%Initializing displacement and force vectors
u = zeros(nDofs,1);
iF = reshape(dim*repmat(loads(:,1),1,dim)-((dim-1):-1:0),[],1);
jF = ones(size(loads,1)*dim,1);
kF = reshape(loads(:,2:end),[],1);
f = sparse(iF,jF,kF,nDofs,1);

nSupp = size(supp,1);
prscDofs = bsxfun(@minus,dim*supp(1:nSupp,1),(dim-1:-1:0)).*supp(1:nSupp,2:end);
prscDofs = prscDofs(prscDofs>0);
freeDofs=setdiff((1:nDofs)',prscDofs);

%Defining filter for sensitivity
rad = filterRad*(dx+dy)/2;
midpoint=[mean(reshape(node(elem',1),4,nElem),1)' mean(reshape(node(elem',2),4,nElem),1)'];
dist1 = sqrt((midpoint(:,1)-midpoint(:,1)').^2+(midpoint(:,2)-midpoint(:,2)').^2);
ins = le(dist1,rad);
[iIn,jIn] = find(ins);      kIn = rad-dist1(ins(:));
H = sparse(iIn,jIn,kIn,nElem,nElem);
Hs = sum(H,2);

%Initializing density and fetch lattice interpolation data
density=volFrac*ones(nElem,1);

ref_EVec = getRefVec(msType,matProp(3));  %Initially extracting reference eigVec
%  Get coeff for chosen degree and lattice
[EigValCoeff,derEigValCoeff] = getEigValCoeff(msType,deg,matProp(2),matProp(3));
[EigAngCoeff,derEigAngCoeff] = getEigAngCoeff(msType);

%  Start of analysis
[Cs,gradCs] = getCinterp(density,EigValCoeff,EigAngCoeff,ref_EVec,sizeOfC,derEigValCoeff,derEigAngCoeff,map);
Ke = getKElem(Ke_template,Cs);
K = sparse(iK,jK,Ke(:),nDofs,nDofs);
F=f-K(:,prscDofs)*u(prscDofs,1);
u(freeDofs,1) = K(freeDofs,freeDofs)\F(freeDofs);
gradKe = getKElem(Ke_template,gradCs);
switch objType
    case 1
        J = u'*f;       %Computing objective and gradients for compliance
        grad = -sum(squeeze(sum(u(elemDofs).*permute(gradKe,[3,1,2]),2)).*u(elemDofs),2);
end

%Initializing storage variables and objective scaling
objHist = zeros(maxIter,1);     change = zeros(maxIter,1);
scaleObj = J;    disp(['Scaling objective by a factor of ' num2str(scaleObj)]);
for iter = 1:maxIter                    %Start iterations
    objHist(iter) = J;             %Scaling objective and gradients
    grad = grad/scaleObj;
    grad_regularized = (H*(density.*grad))./Hs./density;    %Applying sensitivity filter
    %Updating design variable - the densities
    [density,change(iter)] = updateRho1(density,grad_regularized,volFrac*nElem,0.2,1e-6);
    %Displaying result in command
    fprintf('Iter. = %d, obj = %f, change = %f, vol = %f\n',iter,objHist(iter),change(iter),sum(density)*dx*dy);
    if iter>2
        if((abs((objHist(iter)-objHist(iter -1))/objHist(iter-1))<1e-6) && (abs((objHist(iter-1)-objHist(iter -2))/objHist(iter-1))<1e-6))
            disp('Objective Function Converged');   
            break;
        end
    end
    if(change(iter)<0.001)
        disp('Change in density converged');        
        break;
    end
    %Computing interpolated D and its gradient
    [Cs,gradCs] = getCinterp(density,EigValCoeff,EigAngCoeff,ref_EVec,sizeOfC,derEigValCoeff,derEigAngCoeff,map);
    Ke = getKElem(Ke_template,Cs);      %Computing elemental stiffness matrices using templates and interpolants
    K = sparse(iK,jK,Ke(:),nDofs,nDofs);
    F=f-K(:,prscDofs)*u(prscDofs,1);
    u(freeDofs,1) = K(freeDofs,freeDofs)\F(freeDofs);
    gradKe = getKElem(Ke_template,gradCs);
    switch objType
        case 1          %Computing objective and gradients for compliance
            J = u'*f;
            grad = -sum(squeeze(sum(u(elemDofs).*permute(gradKe,[3,1,2]),2)).*u(elemDofs),2);
    end
end
%--------------------------------------------%
function K_elem = getKElem(Ke_template,Ds)
K_elem = zeros(8,8,size(Ds,3));
for comp = 1:size(Ds,1)
    K_elem = K_elem+Ke_template(:,:,comp).*Ds(comp,1,:);
end
%--------------------------------------------%
function [xnew,change] = updateRho1(x,dc,totVol,move,xmin)
%% OPTIMALITY CRITERIA UPDATE OF DESIGN VARIABLES AND PHYSICAL DENSITIES
%Direct Lagrange Multiplier Updates (Kumar & Suresh (2020))
change1 = true;eta = 0.5;ocIter = 0;
nElems = numel(x);dv = ones(nElems,1)*totVol/nElems;
varIn = true(nElems,1);
xMax = min(x+move,1);xMin = max(x-move,xmin);
gRem = totVol; gToDist = gRem;
dc(:) = min(dc(:),0);
xTimesGrad = x.*((-dc./dv).^eta);
while(change1)
    xnew = xTimesGrad/(sum(xTimesGrad(varIn))/gToDist);
    upLgc = xnew>xMax;dnLgc = xnew<xMin;
    gToDist = gRem - sum(xMax(upLgc))-sum(xMin(dnLgc));
    change1 = ~isequal(~varIn,(upLgc|dnLgc));
    varIn = ~(upLgc|dnLgc);ocIter = ocIter + 1;
end
xnew(upLgc) = xMax(upLgc);xnew(dnLgc) = xMin(dnLgc);
change = max(abs(xnew(:)-x(:)));