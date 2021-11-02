function [C_comp] = getCcomp(rho_vol,eigValCoeff,sizeOfC,eVec,map)
nElem = numel(rho_vol);
eVal = zeros(sizeOfC,nElem); %eigval for all elements
for i = 1:sizeOfC
    eVal(i,:) = polyval(eigValCoeff(i,:),rho_vol)';
end
eVal = reshape(eVal,1,sizeOfC,[]);
vL = (eVec.*eVal);%eigvec*eigval: each col of eVec is multiplied with respective eigVal
C =  [sum(vL.*eVec(1,:,:),2),sum(vL.*eVec(2,:,:),2),sum(vL.*eVec(3,:,:),2)];%D(:,i) = vL.*repmat(eVec(i,:),3,1)
C_comp = zeros(sizeOfC*(sizeOfC+1)/2,1,nElem);
for i = 1:sizeOfC*(sizeOfC+1)/2
    C_comp(i,1,:) = C(map{i}(1),map{i}(2),:);
end