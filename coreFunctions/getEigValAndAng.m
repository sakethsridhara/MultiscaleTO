function [eValArr,thetaArr,ref_EVec] = getEigValAndAng(D,density,ref_EVec)
EVec = zeros(size(D));
EVal = zeros(size(D));
instances = size(D,3);
for i = instances:-1:2  %i = 1 is a zero elasticity matrix
    [EVec(:,:,i),EVal(:,:,i)] = eig(D(:,:,i),'balance');
    if i == instances
        temp = diag(EVal(:,:,i));
        [~,ind] = sort(temp,'descend');%Index according to mag. of eig at rho= 1
        EVal(:,:,i) = diag(temp(ind));
        EVec(:,:,i) = EVec(:,ind,i);
    elseif i >= 1
        [~,ind] = max(abs(EVec(:,:,i+1)'*EVec(:,:,i)),[],2);
        EVec(:,:,i) = EVec(:,ind,i);
        DiagonolTerm= diag(EVec(:,:,i+1)'*EVec(:,:,i));
        for j = 1:3
            if DiagonolTerm(j) <1e-3
                EVec(:,j,i) = -EVec(:,j,i);
            end
        end
        temp = diag(EVal(:,:,i));
        EVal(:,:,i) = diag(temp(ind));
    end
end
%Store eigenvalues as vectors instead of matrix
eValArr = zeros(instances,3);
eValArr(:,1) = squeeze(EVal(1,1,:));
eValArr(:,2) = squeeze(EVal(2,2,:));
eValArr(:,3) = squeeze(EVal(3,3,:));

%Use givens rotation
alpha = zeros(instances,1);
beta = zeros(instances,1);
gamma = zeros(instances,1);

if nargin<3
    ref_EVec = EVec(:,:,2);
end

for i = 2:instances     %i = 1 is a zero elasticity matrix
    thisEVec = EVec(:,:,i);
    R = thisEVec/ref_EVec;
    alpha(i) = atan2(R(2,1),R(1,1));
    c = cos(alpha(i));s = sin(alpha(i));
    Rz = blkdiag([c -s;s c],1);
    S = Rz.'*R;
    beta(i) = atan2(-S(3,1),S(1,1));
    c = cos(beta(i));s = sin(beta(i));
    Ry = [c 0 s;0 1 0;-s 0 c];
    T = Ry.'*S;
    gamma(i) = atan2(T(3,2),T(2,2));
end
thetaArr = [alpha,beta,gamma];
%Assign i = 1 (density=0)angles according to the slope from next two values
thetaArr(1,:) = thetaArr(2,:) - density(2)*(thetaArr(3,:)-thetaArr(2,:))/(density(3)-density(2));
