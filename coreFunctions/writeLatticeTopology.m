function writeLatticeTopology(lType)
%WRITELATTICETOPOLOGY 
if nargin <1
    lType = 8;
end

%nodes = [0 0;1 0;1 1;0 1;0.5 0;1 0.5;0.5 1;0 0.5;0.5 0.5];
switch lType 
    case 1  %Box
         nodes = [0 0;1 0;1 1;0 1];   
         conn = [1 2;2 3;3 4;4 1];
    case 2  %Plus
        nodes = [0.5 0;1 0.5;0.5 1;0 0.5];
        conn = [1 3;2 4];
    case 3  %inverted-Z
         nodes = [0 0;1 0;1 1;0 1];   
         conn = [1 2;2 4;3 4];
    case 4  %Z
         nodes = [0 0;1 0;1 1;0 1];   
         conn = [1 2;1 3;3 4];
    case 5  %N
         nodes = [0 0;1 0;1 1;0 1];   
         conn = [2 3;2 4;1 4];
    case 6  %inverted-N
         nodes = [0 0;1 0;1 1;0 1];   
         conn = [2 3;1 3;1 4];
    case 7  %X
         nodes = [0 0;1 0;1 1;0 1];   
         conn = [1 3;2 4];
    case 8  %Hexagonal
         s = 1;x = 3*s;y = sqrt(3)*s; 
         nodes = [s,0;2*s,0;s*2.5,y/2;2*s,y;s,y;s/2,y/2;0,y/2;x,y/2];
         conn = [1,2;2,3;3 4;4 5;5 6;6 1;6 7;3 8];
    case 9 %Kagome = isotropic 
         s = 1;x = 3*s;y = sqrt(3)*s; xf1 = x/4;
         nodes = [0,y/2;x/2,0;x,y/2;x/2,y;xf1,0;xf1,y;x-xf1,0;x-xf1,y];
         conn = [1,2;2,3;3 4;4 1;5 6;7 8];
    case 10 %Diamond = orthotropic
         s = 1;x = 3*s;y = sqrt(3)*s;
         nodes = [0,y/2;x/2,0;x,y/2;x/2,y];
         conn = [1,2;2,3;3 4;4 1;2 4]; 
    case 11 %Triangular honeycomb = almost isotropic
         s = 1;x = 3*s;y = sqrt(3)*s;
         nodes = [0,y/2;x/2,0;x,y/2;x/2,y;0 0;x 0;x y;0 y];
         conn = [1,2;2,3;3 4;4 1;2 4;5 8;6 7]; 
    case 12 %Diamond-Hexagonal = orthotropic
         s = 1;x = 2*s;y = sqrt(3)*s;
         nodes = [0,y/2;x/2,0;x,y/2;x/2,y;x/4,y/2;3*x/4,y/2];
         conn = [5,2;2 6;6 4;4 5;1 5;6 3]; 
    case 13 %Diamond-Octagon = isotropic orthotropic
         s = 1/(sqrt(2) + 1);x = (sqrt(2) + 1)*s;y = (sqrt(2) + 1)*s;f = s/sqrt(2);
         nodes = [0,y/2;x/2,0;x,y/2;x/2,y;x/2-f,y/2;x/2,y/2-f;x/2+f,y/2;x/2,y/2+f];
         conn = [5,6;6 7;7 8;8 5;1 5;2 6;3 7;4 8]; 
    case 14 %Box+X
         nodes = [0 0;1 0;1 1;0 1];   
         conn = [1 2;2 3;3 4;4 1;1 3;2 4];
    otherwise%To write a temporary lattice connectivity        
end
nodes = nodes./min(max(nodes,[],1));   %Ensuring the smaller dimension is unit
folder = '../latticeTopologies/';
filename = ['lattice2d_',num2str(lType)];
ext = '.txt';
address = [folder,filename,ext];
fid = fopen(address,'w');
fprintf(fid,"%i %i\n",max(conn,[],'all'),size(conn,1));
for i = 1:max(conn,[],'all')
    fprintf(fid,"%8.5f %8.5f\n",nodes(i,1),nodes(i,2));
end
for i = 1:size(conn,1)
    fprintf(fid,"%i %i\n",conn(i,1),conn(i,2));
end
fclose(fid);
end