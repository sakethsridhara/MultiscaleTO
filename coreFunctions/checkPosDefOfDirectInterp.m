function posDef = checkPosDefOfDirectInterp(msType,degree,testInstances)
densityTest = linspace(0,1,testInstances+1);
densityTest = densityTest(2:end);
C_int_Direct = zeros(testInstances,6);
folder = './FitData/';filename = 'directInterpCoeff'; suffix = '.mat';
fname = [folder filename num2str(degree*100+msType) suffix];
load(fname,'interpCoeff');
C_temp = zeros(3,3,testInstances);
map = {[1,1],[2,2],[3,3],[1,2],[1,3],[2,3]};
for j = 1:6
    C_int_Direct(:,j) = polyval(interpCoeff(j,:),densityTest);
    C_temp(map{j}(1),map{j}(2),:) = reshape(C_int_Direct(:,j),1,1,[]);
end
posDef = true;
for j = 1:testInstances
    try chol(C_temp(:,:,j));
    catch
        posDef = false;
        break;
    end
end