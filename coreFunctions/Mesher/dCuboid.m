function d = dCuboid(P,x1,x2,y1,y2,z1,z2)
d = [x1-P(:,1), P(:,1)-x2, y1-P(:,2), P(:,2)-y2, z1-P(:,3), P(:,3)-z2];
d = [d,max(d,[],2)];%return dist and max distance