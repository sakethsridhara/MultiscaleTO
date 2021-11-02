function [coeff,derCoeff] = getCdirectCoeff(msType,deg)
filename = 'directInterpCoeff'; suffix = '.mat';
fnameEval = [filename num2str(deg*100+msType)  suffix];
folder = './FitData/';
load([folder fnameEval],'interpCoeff');
coeff = interpCoeff;
if nargout>1
    derCoeff = zeros(size(coeff,1),size(coeff,2)-1);
    for comp = 1:size(coeff,1)
        tempVar = polyder(coeff(comp,:));   %if the higher order coefficients are zero, der may return a lower order polynomial
        derCoeff(comp,:) = [zeros(1,size(coeff,2)-1-length(tempVar)),tempVar];
    end
end