%------------------------------ PolyMesher -------------------------------%
% Ref: C Talischi, GH Paulino, A Pereira, IFM Menezes, "PolyMesher: A     %
%      general-purpose mesh generator for polygonal elements written in   %
%      Matlab," Struct Multidisc Optim, DOI 10.1007/s00158-011-0706-z     %
%-------------------------------------------------------------------------%
function [x] = MbbDomain(Demand,dim,Arg)
switch dim
    case 2
        BdBox = [0 2 0 1];
    case 3
        BdBox = [0 2 0 1 -0.2 0.2];
end
switch(Demand)
    case('Dist');  x = DistFnc(dim,Arg,BdBox);
    case('BC');    x = BndryCnds(dim,Arg{:},BdBox);
    case('BdBox'); x = BdBox;
    case('PFix');  x = FixedPoints(BdBox);
end
%----------------------------------------------- COMPUTE DISTANCE FUNCTIONS
function Dist = DistFnc(dim,P,BdBox)
switch dim
    case 2
        Dist = dRectangle(P,BdBox(1),BdBox(2),BdBox(3),BdBox(4));
    case 3
        Dist = dCuboid(P,BdBox(1),BdBox(2),BdBox(3),BdBox(4),BdBox(5),BdBox(6));
end
%---------------------------------------------- SPECIFY BOUNDARY CONDITIONS
function [x] = BndryCnds(dim,Node,~,BdBox)
switch dim
    case 2
        eps = 0.1*sqrt((BdBox(2)-BdBox(1))*(BdBox(4)-BdBox(3))/size(Node,1));
        LeftEdgeNodes = find(abs(Node(:,1)-BdBox(1))<eps);
        LeftUpperNode = find(abs(Node(:,1)-BdBox(1))<eps & ...
            abs(Node(:,2)-BdBox(4))<eps);
        RigthBottomNode = find(abs(Node(:,1)-BdBox(2))<eps & ...
            abs(Node(:,2)-BdBox(3))<eps);
        FixedNodes = [LeftEdgeNodes; RigthBottomNode];
        Supp = zeros(length(FixedNodes),3);
        Supp(:,1)=FixedNodes; Supp(1:end-1,2)=1; Supp(end,3)=1;
        Load = [LeftUpperNode,0,-1];
    case 3
        eps = 0.1*sqrt((BdBox(2)-BdBox(1))*(BdBox(4)-BdBox(3))*(BdBox(6)-BdBox(5))/size(Node,1));
        LeftEdgeNodes = find(abs(Node(:,1)-BdBox(1))<eps);
        LeftUpperNode = find(abs(Node(:,1)-BdBox(1))<eps & ...
            abs(Node(:,2)-BdBox(4))<eps);
        RigthBottomFrontNode = find(abs(Node(:,1)-BdBox(2))<eps & ...
            abs(Node(:,2)-BdBox(3))<eps & abs(Node(:,3)-BdBox(5))<eps);
        RigthBottomBackNode = find(abs(Node(:,1)-BdBox(2))<eps & ...
            abs(Node(:,2)-BdBox(3))<eps & abs(Node(:,3)-BdBox(6))<eps);
        Supp1 = zeros(length(LeftEdgeNodes),4);Supp1(:,1) = LeftEdgeNodes;Supp1(:,2) = 1;
        Supp2 = zeros(length(RigthBottomFrontNode),4);Supp2(:,1) = RigthBottomFrontNode;        
        Supp3 = zeros(length(RigthBottomBackNode),4);Supp3(:,1) = RigthBottomBackNode;
        Supp2(:,3) = 1;Supp2(:,4) = 1; Supp3(:,3) = 1; 
        Supp=[Supp1;Supp2;Supp3];
        Load = [LeftUpperNode,ones(size(LeftUpperNode))*[0,-0.5,0]];
end
x = {Supp,Load};
%----------------------------------------------------- SPECIFY FIXED POINTS
function [PFix] = FixedPoints(BdBox)
PFix = [];
%-------------------------------------------------------------------------%