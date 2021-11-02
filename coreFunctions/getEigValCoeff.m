function [coeff,derCoeff] = getEigValCoeff(Lattice_type,degreeOfInterp,YM,nu)
if nargin<3,    YM = 1; end
if nargin<4,    nu = 0.3; end
if Lattice_type == 2
    latToGet = 1;
else
    latToGet = Lattice_type;
end
lattice = degreeOfInterp*100+latToGet; %identifier for filename
if Lattice_type == 0 %SIMP is MS_type = 0
    if nargin<5,        nu = 0.3;    end
    D = 1/(1-nu^2)*[1 nu 0;nu 1 0;0 0 (1-nu)/2];
    eVals = eig(D);
    evalCoeff = eVals.*[1,0,0,0];
    rotCoeff = zeros(3,size(evalCoeff,2));
    coeff = [evalCoeff;rotCoeff];
else
    filename = 'eigValCoeff'; suffix = '.mat';
    fnameEval = [filename num2str(lattice) suffix];
    folder = './FitData/';
    load([folder fnameEval],'eValCoeff');
    coeff = eValCoeff*YM;
end
derCoeff = zeros(size(coeff,1),size(coeff,2)-1);
for comp = 1:size(coeff,1)
    tempVar = polyder(coeff(comp,:));   %if the higher order coefficients are zero, der may return a lower order polynomial
    derCoeff(comp,:) = [zeros(1,size(coeff,2)-1-length(tempVar)),tempVar];
end