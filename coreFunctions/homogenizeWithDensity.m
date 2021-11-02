function CH = homogenizeWithDensity(lx, ly, lambda, mu, phi, x, p)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% lx        = Unit cell length in x-direction.
% ly        = Unit cell length in y-direction.
% lambda    = Lame's first parameter for both materials. Two entries.
% mu        = Lame's second parameter for both materials. Two entries.
% phi       = Angle between horizontal and vertical cell wall. Degrees
% x         = Material (2) fraction matrix. Size used to determine nelx/nely
% p         = Penalty        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INITIALIZE
% Deduce discretization: x 
[nely, nelx] = size(x);
% Stiffness matrix consists of two parts, one belonging to lambda and
% one belonging to mu. Same goes for load vector
dx = lx/nelx; dy = ly/nely;
nel = nelx*nely;
[keLambda, keMu, feLambda, feMu] = elementMatVec(dx/2, dy/2, phi);
% Node numbers and element degrees of freedom for full (not periodic) mesh
nodenrs = reshape(1:(1+nelx)*(1+nely),1+nely,1+nelx);
edofVec = reshape(2*nodenrs(1:end-1,1:end-1)+1,nel,1);
edofMat = repmat(edofVec,1,8)+repmat([0 1 2*nely+[2 3 0 1] -2 -1],nel,1);
%% IMPOSE PERIODIC BOUNDARY CONDITIONS
% Use original edofMat to index into list with the periodic dofs
nn = (nelx+1)*(nely+1); % Total number of nodes
nnP = (nelx)*(nely);    % Total number of unique nodes
nnPArray = reshape(1:nnP, nely, nelx);
% Extend with a mirror of the top border
nnPArray(end+1,:) = nnPArray(1,:);
% Extend with a mirror of the left border
nnPArray(:,end+1) = nnPArray(:,1);
% Make a vector into which we can index using edofMat:
dofVector = zeros(2*nn, 1);
dofVector(1:2:end) = 2*nnPArray(:)-1;
dofVector(2:2:end) = 2*nnPArray(:);
edofMat = dofVector(edofMat);
ndof = 2*nnP; % Number of dofs
%% ASSEMBLE STIFFNESS MATRIX
% Indexing vectors
iK = kron(edofMat,ones(8,1))';
jK = kron(edofMat,ones(1,8))';
% Material properties in the different elements
lambdaS = lambda(1)*(1-x.^p) + lambda(2)*(x.^p);
muS     = mu(1)*(1-x.^p) + mu(2)*(x.^p);
% The corresponding stiffness matrix entries
sK = keLambda(:)*lambdaS(:).' + keMu(:)*muS(:).';
K  = sparse(iK(:), jK(:), sK(:), ndof, ndof); 
%% LOAD VECTORS AND SOLUTION
% Assembly three load cases corresponding to the three strain cases
sF = feLambda(:)*lambdaS(:).'+feMu(:)*muS(:).';
iF = repmat(edofMat',3,1);
jF = [ones(8,nel); 2*ones(8,nel); 3*ones(8,nel)];
F  = sparse(iF(:), jF(:), sF(:), ndof, 3);
% Solve (remember to constrain one node)
chi(3:ndof,:) = K(3:ndof,3:ndof)\F(3:ndof,:);
%% HOMOGENIZATION
% The displacement vectors corresponding to the unit strain cases
chi0 = zeros(nel, 8, 3);
% The element displacements for the three unit strains
chi0_e = zeros(8, 3);
ke = keMu + keLambda; % Here the exact ratio does not matter, because
fe = feMu + feLambda; % it is reflected in the load vector
chi0_e([3 5:end],:) = ke([3 5:end],[3 5:end])\fe([3 5:end],:);
% epsilon0_11 = (1, 0, 0)
chi0(:,:,1) = kron(chi0_e(:,1)', ones(nel,1));
% epsilon0_22 = (0, 1, 0)
chi0(:,:,2) = kron(chi0_e(:,2)', ones(nel,1));
% epsilon0_12 = (0, 0, 1)
chi0(:,:,3) = kron(chi0_e(:,3)', ones(nel,1));
CH = zeros(3);
cellVolume = lx*ly;
for i = 1:3
  for j = 1:3
    sumLambda = ((chi0(:,:,i) - chi(edofMat+(i-1)*ndof))*keLambda).*...
      (chi0(:,:,j) - chi(edofMat+(j-1)*ndof));
    sumMu = ((chi0(:,:,i) - chi(edofMat+(i-1)*ndof))*keMu).*...
      (chi0(:,:,j) - chi(edofMat+(j-1)*ndof));
    sumLambda = reshape(sum(sumLambda,2), nely, nelx);
    sumMu = reshape(sum(sumMu,2), nely, nelx);
    % Homogenized elasticity tensor
    CH(i,j) = 1/cellVolume*sum(sum(lambdaS.*sumLambda + muS.*sumMu));
  end
end
%% COMPUTE ELEMENT STIFFNESS MATRIX AND FORCE VECTOR (NUMERICALLY)
function [keLambda, keMu, feLambda, feMu] = elementMatVec(a, b, phi)
% Constitutive matrix contributions
CMu = diag([2 2 1]); CLambda = zeros(3); CLambda(1:2,1:2) = 1; 
% Two Gauss points in both directions
xx=[-1/sqrt(3), 1/sqrt(3)]; yy = xx;
ww=[1,1];
% Initialize
keLambda = zeros(8,8); keMu = zeros(8,8);
feLambda = zeros(8,3); feMu = zeros(8,3);
L = zeros(3,4); L(1,1) = 1; L(2,4) = 1; L(3,2:3) = 1;
for ii=1:length(xx)
  for jj=1:length(yy)
    % Integration point
    x = xx(ii); y = yy(jj);
    % Differentiated shape functions
    dNx = 1/4*[-(1-y)  (1-y) (1+y) -(1+y)];
    dNy = 1/4*[-(1-x) -(1+x) (1+x)  (1-x)];
    % Jacobian
    J = [dNx; dNy]*[-a a a+2*b/tan(phi*pi/180) 2*b/tan(phi*pi/180)-a; ...
        -b -b b b]';
    detJ = J(1,1)*J(2,2) - J(1,2)*J(2,1);
    invJ = 1/detJ*[J(2,2) -J(1,2); -J(2,1) J(1,1)];
    % Weight factor at this point
    weight = ww(ii)*ww(jj)*detJ;
    % Strain-displacement matrix
    G = [invJ zeros(2); zeros(2) invJ];
    dN = zeros(4,8);
    dN(1,1:2:8) = dNx;
    dN(2,1:2:8) = dNy;
    dN(3,2:2:8) = dNx;
    dN(4,2:2:8) = dNy;
    B = L*G*dN;
    % Element matrices
    keLambda = keLambda + weight*(B' * CLambda * B);
    keMu = keMu + weight*(B' * CMu * B);
    % Element loads
    feLambda = feLambda + weight*(B' * CLambda * diag([1 1 1]));
    feMu = feMu + weight*(B' * CMu * diag([1 1 1]));
  end
end