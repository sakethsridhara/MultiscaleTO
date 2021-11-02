%% Sec1: Run Desired Sections use either "Run Section" or "Run and Advance"
clear; clc; close all;
%% Sec2: Adding all subfolders to path
fList = dir;flag = false;
for i = 3:size(fList,1)
    if fList(i).isdir && strcmp('FitData',fList(i).name)
        flag = true;
        break;
    end
end
if flag == false
    mkdir FitData;
end
folder = fileparts(which(mfilename));
addpath(genpath(folder));
%% Sec3: Set parameters here
% List of microstructure types
% 1 = square, 2 = cross, 3 = inverted-Z; 4 = Z; 5 = N; 6 = inverted-N; 
% 7 = X; 8 = hexagonal; 9 = kagome; 10 = mixed-triangular; 
% 11 = triangular-hexagonal; 12 = diamond-hexagonal
msType = 7;
degree = 4;         %degree of polynomial fitting function

%Optimization domain: 1 = Cantilever, 2 = LBracket, 3 = MBB
domainType = 3;
volFrac = 0.5;

%% Sec4: Define default parameters
nely = 100;                 %Microstructure discretization; nelx = AR*nely
nu = 0.3;                   %Base material Poisson's ratio nu; E = 1
instances = 11;             %Instances = s+1; s>=degree


%% Sec5: CurveFitting for Spectral Decomposition (SD)
disp('Computing fitting functions ...');
generateFitFuncsForSD(msType,nu,nely,instances,degree);
%% Sec6: TopOpt using Spectral Decomposition (SD)
OptimizeSD(domainType, msType, degree, volFrac);
%% Sec7: CurveFitting for Direct Regression (SD) along with pos-def check
disp('Computing fitting functions ...')
generateFitFuncsForDirect(msType,nu,nely,instances,degree);
isPosDef = checkPosDefOfDirectInterp(msType,degree,1000);
%% Sec8: TopOpt using Direct Regression (SD)
if isPosDef
    OptimizeDirect(domainType, msType, degree, volFrac);
else
    disp('Non-positive-definite matrix detected; Optimization not possible');
end
%% Sec9: Compare R2 for the two methods
testInstances = 50;
compareR2(msType,nu,nely,testInstances,degree);