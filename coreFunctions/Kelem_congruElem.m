function [K] = Kelem_congruElem(D,t,B,Det)
wgt = [1.0  1.0  1.0  1.0];
K=zeros(8,8);
for i=1:4
    wt = wgt(i)*t*Det;
    K = K + (B{1,i})'*D*B{1,i}*wt;
end