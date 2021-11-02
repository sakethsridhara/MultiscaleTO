function AR = getAspectRatio(msType)
folder = './latticeTopologies/';
filename = ['lattice2d_',num2str(msType)];
ext = '.txt';
address = [folder,filename,ext];
[nodeList] = ReadTrussTopology(address,2);
AR = (max(nodeList(:,1))-min(nodeList(:,1)))/(max(nodeList(:,2))-min(nodeList(:,2)));
end