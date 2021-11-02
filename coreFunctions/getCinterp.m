function [Cs,gradCs] = getCinterp(rho_vol,eigValCoeff,eigAngCoeff,ref_EVec,sizeOfD,derEigValCoeff,derEigAngCoeff,map)
[eVec,derEVec] = getRotEVec(rho_vol,eigAngCoeff,ref_EVec,sizeOfD,derEigAngCoeff);
Cs = getCcomp(rho_vol,eigValCoeff,sizeOfD,eVec,map);
if nargout>1
    gradCs = getDerCcomp(rho_vol,eigValCoeff,derEigValCoeff,sizeOfD,eVec,map,derEVec);
end