function [density,xs,ds] = getDensityAndDistribution(instances,MSType,nely,saveThickness)
%% Find density distribution and compute density
folder = './latticeTopologies/';ext = '.txt';

filename = ['lattice2d_',num2str(MSType)];
address = [folder,filename,ext];
AR = getAspectRatio(MSType);
nels = [nely*AR,nely];

%Find minimium possible density
[~,~,ds] = GeneratePixel5(nels,address,1);
min_t = 2*max(ds./nels);
[~,min_rho] = GeneratePixel5(nels,address,min_t);

%Allocate memory for a 0 density as well
desiredRho = [0;linspace(min_rho,1,instances-1)'];
density = zeros(size(desiredRho));
thickness = zeros(size(desiredRho));
flag = 0;
for i = 1:length(desiredRho)
    if desiredRho(i) == 0, continue; end
    t1 = min_t;t2 = 1;
    while (t2-t1)>1e-6
        tmid = (t2+t1)/2;
        [x,mi_rho] = GeneratePixel5(nels,address,tmid);
        if mi_rho >= desiredRho(i), t2 = tmid; else, t1 = tmid; end
    end
    if flag == 0
        xs = zeros([size(x),length(desiredRho)]);
        flag = 1;
    end
    thickness(i) = tmid;
    xs(:,:,i) = x;
    density(i) = mi_rho;
end
%% Store thickness vs density data to aid in printing optimized designs
if nargin>3
    if saveThickness
        interpType = @pchip;
        folder2 = './FitData/';
        fileToWrite = 'thickVsDensity';
        suffix = '.mat';
        tfunc = interpType(density,thickness);
        fname = [folder2 fileToWrite num2str(MSType) suffix];
        save(fname,'tfunc');
    end
end

