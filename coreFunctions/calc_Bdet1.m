function [B,Det]=calc_Bdet1(Elem,Node,isCongruent)
if nargin <3
    isCongruent = 0;
end
nElem=size(Elem,1);
B=cell(nElem,4);
Det=cell(nElem,4);
XG = [-1/sqrt(3) 1/sqrt(3)  1/sqrt(3) -1/sqrt(3) ];
YG = [-1/sqrt(3) -1/sqrt(3) 1/sqrt(3)  1/sqrt(3)];
XX=zeros(4,2);
for i=1:nElem
    XX(:,1) =  Node(Elem(i,:),1); XX(:,2) =  Node(Elem(i,:),2);
    for GSP = 1:4
        RI = XG(GSP);        SI = YG(GSP);
        [B1, Det1] = B_Det(XX,RI,SI);
        B{i,GSP}=B1;
        Det{i,GSP}=Det1;
    end
    if isCongruent
        B = repmat(B(1,:),nElem,1);
        Det = repmat(Det(1,:),nElem,1);
        break;
    end
end
%--------------------------------------------%
function [B0, Det] = B_Det(XX,R,S)
ZP = 1.0+R;     EP = 1.0+S;   % zeta plus 1 %eta puls 1
ZM = 1.0-R;     EM = 1.0-S;
P = 0.25*[-EM EM EP -EP;-ZM -ZP ZP ZM];
%Evaluation of the Jacobian matrix
XJ  =P*XX;
% Determinant
Det = XJ(1,1)*XJ(2,2)-XJ(1,2)*XJ(2,1);
if(Det <= 0.0)
    disp('ERROR: Jacobian matrix is not positive definite.');
    return;
end
%inverse of jacobian matrix
XJI=zeros(2,2);
XJI(1,1) = XJ(2,2)/Det;
XJI(1,2) = -XJ(1,2)/Det;
XJI(2,1) = -XJ(2,1)/Det;
XJI(2,2) = XJ(1,1)/Det;

R1=[1 0 0 0;0 0 0 1;0 1 1 0];
R2=zeros(4);    R2(1:2,1:2)=XJI;    R2(3:4,3:4)=XJI;
R3=zeros(4,8);  R3(1:2,1:2:7)=P;    R3(3:4,2:2:8)=P;
B0=R1*R2*R3;