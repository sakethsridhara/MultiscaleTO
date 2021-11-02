function compareR2(msType,nu,nely,testInstances,degree)
%% For test instances get density
%[densityTest] = linspace(0,1,testInstances)';
%[densityTest] = getDensityAndDistribution(testInstances,msType,nelx,nely,0);
%% Compute C using AH
[~,density,C_AH_vec] = getHomogenizedData(testInstances,nu,msType,nely);
C_AH_vec = C_AH_vec';
%% Get C from SD fitting functions
sizeOfC = 3;
ref_EVec = getRefVec(msType);  %Initially extracting reference eigVec

% Get coefficient for chosen degree and lattice
[EigValCoeff,derEigValCoeff] = getEigValCoeff(msType,degree);
[EigAngCoeff,derEigAngCoeff] = getEigAngCoeff(msType);

map = {[1,1],[2,2],[3,3],[1,2],[1,3],[2,3]};
[C_int_SD] = getCinterp(density,EigValCoeff,EigAngCoeff,ref_EVec,sizeOfC,derEigValCoeff,derEigAngCoeff,map);
C_int_SD = permute(C_int_SD,[1,3,2]);

%% Now carryout R2 value computation for SD method
R2_SD = zeros(6,1);
for k = 1:6
    R2_SD(k) = rsquare(C_AH_vec(k,:),C_int_SD(k,:));
end
disp('R^2 for SD method');
disp(num2str(R2_SD'));

%% Get C from Direct fitting functions
coeff = getCdirectCoeff(msType,degree);
C_int_Direct = getCdirectInterp(density,coeff);
C_int_Direct = permute(C_int_Direct,[1,3,2]);
%% Now carryout R2 value computation for Direct method
R2_Direct = zeros(6,1);
for k = 1:6
    R2_Direct(k) = rsquare(C_AH_vec(k,:),C_int_Direct(k,:));
end
disp('R^2 for Direct method');
disp(num2str(R2_Direct'));