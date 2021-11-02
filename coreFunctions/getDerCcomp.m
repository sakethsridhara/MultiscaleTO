function [derC_comp] = getDerCcomp(rho_vol,eigValCoeff,derEigValCoeff,sizeOfC,eVec,map,derEVec)
nElem = numel(rho_vol);
eVal = zeros(sizeOfC,nElem);
derEVal = zeros(sizeOfC,nElem);
for i = 1:sizeOfC
    eVal(i,:) = polyval(eigValCoeff(i,:),rho_vol)';
    derEVal(i,:) = polyval(derEigValCoeff(i,:),rho_vol)';
end

derEVal = reshape(derEVal,1,sizeOfC,[]);
vL = (eVec.*derEVal);%each col of eVec is multiplied with respective derivative of eigVal
%Computing second term of derDH = v times derEVal times v transpose
v_derEVal_vt =  [sum(vL.*eVec(1,:,:),2),sum(vL.*eVec(2,:,:),2),sum(vL.*eVec(3,:,:),2)];%matrix(:,i) = vL.*repmat(eVec(i,:),3,1)

eVal = reshape(eVal,1,sizeOfC,[]);
der_vL = (derEVec.*eVal);%der_eVec*val
derv_eVal_vt =  [sum(der_vL.*eVec(1,:,:),2),sum(der_vL.*eVec(2,:,:),2),sum(der_vL.*eVec(3,:,:),2)];
derv_eVal_vt = derv_eVal_vt+permute(derv_eVal_vt,[2,1,3]); %Adding first and third term(transpose of first)
derC = v_derEVal_vt + derv_eVal_vt;

derC_comp = zeros(sizeOfC*(sizeOfC+1)/2,1,nElem);
for i = 1:sizeOfC*(sizeOfC+1)/2
    derC_comp(i,1,:) = derC(map{i}(1),map{i}(2),:);
end