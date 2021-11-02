function [vol,xs,specs] = getDensityAndSpecsConformal(densities,lx,ly,nels,phi)
instances = numel(densities);
vol = zeros(instances,1);
xs = zeros([nels([2,1]),instances]);
specs = zeros(2,4,instances);
ratio_wToh = 1/sind(phi); %phi cannot be zero; 45<=phi<=90
%Assuming lx = ly
w = (ratio_wToh/2)*((1+1/ratio_wToh)-sqrt((1+1/ratio_wToh)^2-4*densities/ratio_wToh));
ws(:,1) = w*(lx+ly)/2;
ws(:,2) = (w/ratio_wToh)*(lx+ly)/2;


nelx = nels(1);nely = nels(2);
for i = 1:instances
    nelby = max(1,floor(nely*min(ly,ws(i,2))/ly));nelty = nely-nelby;
    nellx = max(1,floor(nelx*min(lx,ws(i,1))/lx));nelrx = nelx-nellx;
    
    x = 1*ones(nely,nelx);x(1:nelby,:) = 2;x(:,1:nellx) = 2;
    
    dlx = min(lx,ws(i,1))/nellx; dby = min(ly,ws(i,2))/nelby;
    drx = (lx-min(lx,ws(i,1)))/nelrx; dty = (ly-min(ly,ws(i,2)))/nelty;
    if ~isfinite(dlx), dlx = 0; end;if ~isfinite(dby), dby = 0; end
    if ~isfinite(drx), drx = 0; end;if ~isfinite(dty), dty = 0; end
    xs(:,:,i) = x;
    specs(:,:,i) = [dlx,drx,dby,dty;nellx,nelrx,nelby,nelty];
    vol(i) = nelby*dby*lx+nellx*nelty*dlx*dty+1e-3*nelty*nelrx*drx*dty;
end