function [coeff,derCoeff] = getEigAngCoeff(Lattice_type)
if Lattice_type == 2
    latToGet = 1;
else
    latToGet = Lattice_type;
end
lattice = latToGet;
if Lattice_type == 0 %SIMP is MS_type = 0
    rotCoeff = zeros(3,3);
    coeff = mat2cell(rotCoeff,ones(3,1),3);
else
    filename = 'eigAngCoeff'; suffix = '.mat';
    fnameEval = [filename num2str(lattice) suffix];
    folder = './FitData/';
    load([folder fnameEval],'thetaFit');
    coeff = thetaFit;
end
derCoeff = cell(size(coeff));
for comp = 1:size(coeff,1)
    derCoeff{comp,1} = fnder(coeff{comp,1});
end