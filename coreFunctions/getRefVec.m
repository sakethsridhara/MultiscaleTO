function ref_EVec = getRefVec(Lattice_type,nu)
if Lattice_type == 2
    latToGet = 1;
else
    latToGet = Lattice_type;
end
if Lattice_type == 0 %SIMP is MS_type = 0
    if nargin<2, nu = 0.3; end
    D = 1/(1-nu^2)*[1 nu 0;nu 1 0;0 0 (1-nu)/2];
    [eVecs,~] = eig(D);
    ref_EVec = eVecs;
else
    filename = 'ref_Evec'; suffix = '.mat';
    fnameEvec = [filename num2str(latToGet) suffix];
    folder = './FitData/';
    load([folder fnameEvec],'ref_EVec');%Change name of the variable
end