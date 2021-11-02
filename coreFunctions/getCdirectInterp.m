function [C_comp,derC_comp] = getCdirectInterp(density,msType,deg)
if nargin == 3&&nargout==1
    filename = 'directInterpCoeff'; suffix = '.mat';
    fnameEval = [filename num2str(deg*100+msType)  suffix];
    folder = './FitData/';
    load([folder fnameEval],'interpCoeff');
    coeff = interpCoeff;
elseif nargout == 2
    coeff = msType;
    derCoeff = deg;
elseif nargout == 1 && nargin == 2
    coeff = msType;    
end
nElem = numel(density);
C_comp = zeros(size(coeff,1),nElem); %eigval for all elements
for i = 1:size(coeff,1)
    C_comp(i,:) = polyval(coeff(i,:),density)';
end
C_comp = reshape(C_comp,size(coeff,1),1,nElem);
if nargout>1
    derC_comp = zeros(size(derCoeff,1),nElem); %eigval for all elements
    for i = 1:size(derCoeff,1)
        derC_comp(i,:) = polyval(derCoeff(i,:),density)';
    end
    derC_comp = reshape(derC_comp,size(derCoeff,1),1,nElem);    
end