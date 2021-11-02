function [C_AH,density,C_AH_vec] = getHomogenizedData(instances,nu,MSType,nely,saveThickness)

lambda = nu/((1+nu)*(1-2*nu));mu = 1/(2*(1+nu));
lambdaP = 2*mu*lambda/(lambda+2*mu);
epsilon = 1e-9;
twoMat = [epsilon,1+epsilon];
penalty = 1;

if nargin<5
    saveThickness = 0;
end

%Get uniform density values
[density,xs,ds] = getDensityAndDistribution(instances,MSType,nely,saveThickness);
C_AH = zeros(3,3,length(density));
for i = 1:length(density)
    if density(i) == 0, continue; end
    x = xs(:,:,i);
    [C_AH(:,:,i)] = homogenizeWithDensity(ds(1),ds(2),lambdaP*twoMat,mu*twoMat,90,x,penalty);
end
%get C for rho = 0;
x(:) = 0;
C_eps = homogenizeWithDensity(ds(1),ds(2),lambdaP*twoMat,mu*twoMat,90,x,penalty);

%Add zero density C to data
zeroDensity = density == 0; if sum(zeroDensity)>1, disp('error'); end
C_AH(:,:,zeroDensity) = repmat(C_eps,1,1,sum(zeroDensity));

%Ensure C(0) = 0;
C_AH = C_AH-C_eps;

if nargout>2
    C_AH_vec = zeros(length(density),6);
    C_AH_vec(:,1) = squeeze(C_AH(1,1,:));
    C_AH_vec(:,2) = squeeze(C_AH(2,2,:));
    C_AH_vec(:,3) = squeeze(C_AH(3,3,:));
    C_AH_vec(:,4) = squeeze(C_AH(1,2,:));
    C_AH_vec(:,5) = squeeze(C_AH(1,3,:));
    C_AH_vec(:,6) = squeeze(C_AH(2,3,:));
end
