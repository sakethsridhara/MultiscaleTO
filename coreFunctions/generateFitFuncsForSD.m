function generateFitFuncsForSD(msType,nu,nely,instances,degree)
%% Get Homogenized values
[C_AH,density] = getHomogenizedData(instances,nu,msType,nely,1);
%Calculate eigenvalue and eigenvectors and order them appropriately
[eValArr,thetaArr,ref_EVec] = getEigValAndAng(C_AH,density);

%% Interpolate eigenvalues using Vandermonde matrix
num_eigVals = size(eValArr,2);
num_orientAng = size(thetaArr,2);
folder = './FitData/';
%Initialize interpolation variables
eValCoeff = zeros(num_eigVals,degree+1);
VanderMondeMat = (density.^(degree:-1:0));
p = zeros(num_eigVals,degree+1);
optimOpt = optimoptions('lsqlin','Algorithm','trust-region-reflective',...
    'FunctionTolerance',1.0000e-15,'MaxIterations',1000,'Display','off');
for j = 1:num_eigVals
    eValCoeff(j,1:end-1) = lsqlin(VanderMondeMat(1:end,1:end-1),eValArr(1:end,j),[],[],[],[],zeros(degree,1),inf(degree,1),p(j,1:end-1),optimOpt);
end
%% Interpolate orientation angles using spline interpolation
thetaFit = cell(num_orientAng,1);
for j = 1:num_orientAng
    thetaFit{j,1} = spline(density,thetaArr(:,j));
end
%% Save data
%Save eigenvalue
filename = 'eigValCoeff'; suffix = '.mat';
fname = [folder filename num2str(degree*100+msType) suffix];
save(fname,'eValCoeff');  
%Save thetaFit
filename = 'eigAngCoeff'; suffix = '.mat';
fname = [folder filename num2str(msType) suffix];
save(fname,'thetaFit');
%Save reference eigenvector
filename = 'ref_EVec'; suffix = '.mat';
fname = [folder filename num2str(msType) suffix];
save(fname,'ref_EVec');
