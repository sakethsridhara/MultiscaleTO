function [eVec,derEVec] = getRotEVec(rho_vol,eigAngCoeff,ref_EVec,sizeOfD,derEigAngCoeff)
nElem = numel(rho_vol);
ang = zeros(sizeOfD,nElem);
derAng = zeros(sizeOfD,nElem);
if(isempty(eigAngCoeff))
    eVec = ref_EVec;
    derEVec = zeros(size(ref_EVec));
else
    for i = 1:sizeOfD
        ang(i,:) = ppval(eigAngCoeff{i,1},rho_vol)';
        derAng(i,:) = ppval(derEigAngCoeff{i,1},rho_vol)';
    end
    
    ang = reshape(ang,sizeOfD,1,[]);
    ca = cos(ang(1,1,:));sa = sin(ang(1,1,:));
    cb = cos(ang(2,1,:));sb = sin(ang(2,1,:));
    cg = cos(ang(3,1,:));sg = sin(ang(3,1,:));
    rotMat = [ca.*cb, ca.*sb.*sg-sa.*cg, ca.*sb.*cg+sa.*sg;...
        sa.*cb, sa.*sb.*sg+ca.*cg, sa.*sb.*cg-ca.*sg;...
        -sb, cb.*sg, cb.*cg];
    eVec = permute([sum(rotMat(1,:,:).*ref_EVec',2),...
        sum(rotMat(2,:,:).*ref_EVec',2),sum(rotMat(3,:,:).*ref_EVec',2)],[2,1,3]);
    if nargout>1
        derAng = reshape(derAng,sizeOfD,1,[]);
        dAl = derAng(1,1,:);dBe = derAng(2,1,:);dGa = derAng(3,1,:);
        drotMat = [-dAl.*cb.*sa-dBe.*ca.*sb, dGa.*sa.*sg-dAl.*ca.*cg+dBe.*ca.*cb.*sg+dGa.*ca.*cg.*sb-dAl.*sa.*sb.*sg, dAl.*ca.*sg+dGa.*cg.*sa+dBe.*ca.*cb.*cg-dAl.*cg.*sa.*sb-dGa.*ca.*sb.*sg;...
            dAl.*ca.*cb-dBe.*sa.*sb, dAl.*ca.*sb.*sg-dGa.*ca.*sg-dAl.*cg.*sa+dBe.*cb.*sa.*sg+dGa.*cg.*sa.*sb, dAl.*sa.*sg-dGa.*ca.*cg+dAl.*ca.*cg.*sb+dBe.*cb.*cg.*sa-dGa.*sa.*sb.*sg;...
            -dBe.*cb,                                                dGa.*cb.*cg-dBe.*sb.*sg,                                                     -dBe.*cg.*sb-dGa.*cb.*sg];
        derEVec = permute([sum(drotMat(1,:,:).*ref_EVec',2),...
            sum(drotMat(2,:,:).*ref_EVec',2),sum(drotMat(3,:,:).*ref_EVec',2)],[2,1,3]);
    end
end