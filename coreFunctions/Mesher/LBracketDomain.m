%------------------------------ PolyMesher -------------------------------%
% Ref: C Talischi, GH Paulino, A Pereira, IFM Menezes, "PolyMesher: A     %
%      general-purpose mesh generator for polygonal elements written in   %
%      Matlab," Struct Multidisc Optim, DOI 10.1007/s00158-011-0706-z     %
%-------------------------------------------------------------------------%
function [x] = LBracketDomain(Demand,dim,Arg)
if dim == 2
  BdBox = [0 2.5 0 2.5];
  switch(Demand)
    case('Dist');  x = DistFnc(Arg,BdBox);
    case('BC');    x = BndryCnds(Arg{:},BdBox);
    case('BdBox'); x = BdBox;
    case('PFix');  x = FixedPoints(BdBox);
  end
end
%----------------------------------------------- COMPUTE DISTANCE FUNCTIONS
function Dist = DistFnc(P,BdBox)
  fx = 0.4; fy = 0.4;
  dR1 = dRectangle(P,BdBox(1),BdBox(2),BdBox(3),BdBox(4));
  dR2 = dRectangle(P,BdBox(1)+fx*(BdBox(2)-BdBox(1)),BdBox(2),BdBox(3)+fy*(BdBox(4)-BdBox(3)),BdBox(4));
  Dist = dDiff(dR1,dR2);
%   d1 = dLine(P,0,0.3,0,-0.3);
%   d2 = dLine(P,0,-0.3,2,-0.5);
%   d3 = dLine(P,2,-0.5,2,0.5);
%   d4 = dLine(P,2,0.5,0,0.3);
%   d5 = dCircle(P,0,0,0.3);
%   d6 = dCircle(P,2,0,0.5);
%   douter = dUnion(d6,dUnion(d5,...
%            dIntersect(d4,dIntersect(d3,dIntersect(d2,d1)))));
%   d7 = dCircle(P,0,0,0.175);
%   d8 = dCircle(P,2,0,0.3);
%   din = dUnion(d8,d7);
%   Dist = dDiff(douter,din);
%---------------------------------------------- SPECIFY BOUNDARY CONDITIONS
function [x] = BndryCnds(Node,Element,BdBox)
 % eps = 0.1*sqrt((BdBox(2)-BdBox(1))*(BdBox(4)-BdBox(3))/size(Node,1));
  epsi = 1e-8;
  TopNodes = find(Node(:,2)>max(Node(:,2))-epsi);   
  
  Supp = ones(size(TopNodes,1),3);
  Supp(:,1) = TopNodes;
  RightNodes = find(Node(:,1)>max(Node(:,1))-epsi);   
  [~,TopOfRightNodes] = max(Node(RightNodes,2));
  Load = 1*ones(size(RightNodes(TopOfRightNodes),1),3);
  Load(:,1) = RightNodes(TopOfRightNodes); Load(:,2) = 0;
  x = {Supp,Load};
%----------------------------------------------------- SPECIFY FIXED POINTS
function [PFix] = FixedPoints(BdBox)
  PFix = [];
%-------------------------------------------------------------------------%