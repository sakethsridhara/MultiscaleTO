function generateFitFuncsForDirect(msType,nu,nely,instances,degree)
[~,density,D_vec] = getHomogenizedData(instances,nu,msType,nely,1);
%Interpolate homogenized matrix components directly
num_comp = size(D_vec,2);
folder = './FitData/';

%Initialize interpolation variables
interpCoeff = zeros(num_comp,degree+1);
VanderMondeMat = (density.^(degree:-1:0));
p = zeros(num_comp,degree+1);
optimOpt = optimoptions('lsqlin','Algorithm','trust-region-reflective',...
    'FunctionTolerance',1.0000e-15,'MaxIterations',1000,'Display','off');
for j = 1:num_comp
    interpCoeff(j,1:end-1) = lsqlin(VanderMondeMat(1:end,1:end-1),D_vec(1:end,j),[],[],[],[],[],[],p(j,1:end-1),optimOpt);
end
%Save data to file
filename = 'directInterpCoeff'; suffix = '.mat';
fname = [folder filename num2str(degree*100+msType) suffix];
save(fname,'interpCoeff');


