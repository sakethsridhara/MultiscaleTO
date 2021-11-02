function [nodeList,trussList] = ReadTrussTopology(address,dof)
fid = fopen(address,'r');
nums=fscanf(fid,"%i\n",2);
nodeList = fscanf(fid,"%f %f\n",[dof,nums(1)])';
trussList = fscanf(fid,"%i %i\n",[2,nums(2)])';
fclose(fid);
end