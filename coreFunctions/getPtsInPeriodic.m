function [ptsIn,signedDist] = getPtsInPeriodic(pts,radius,connectivity,node)
ptsIn = false(size(pts,1),1);%ptsIn = inf(size(pts,1),1);
dof = size(pts,2);
signedDist = Inf(size(ptsIn));
tol = 1e-10;
offs = -1:1;
mins = min(node,[],1);maxs = max(node,[],1);

ptsPeriodic = repmat(pts,9,1)...
    +[repelem(repelem([-1,0,1]'*(maxs(1)-mins(1)),3,1),size(pts,1),1),...
    repelem(repmat([-1,0,1]'*(maxs(2)-mins(2)),3,1),size(pts,1),1)];
for j = 1:size(connectivity,1)
    start_n = node(connectivity(j,1),:);  % start node coordinate
    end_n = node(connectivity(j,2),:);    % end node coordinate
    switch dof
        case 2
            % Return minimum distance between line segment vw and point p
            v1 = (end_n-start_n);
            l2 = sum(v1.^2);  % i.e. |w-v|^2 -  avoid a sqrt
            if (l2 == 0.0)
                distance = sum((ptsPeriodic-start_n).^2);   % v == w case
         %       distance = sum((pts-start_n).^2);   % v == w case
            else
                % Consider the line extending the segment, parameterized as v + t (w - v).
                % We find projection of point p onto the line.
                % It falls where t = [(p-v) . (w-v)] / |w-v|^2
                % We clamp t from [0,1] to handle points outside the segment vw.
              
                t = max(0, min(1, (ptsPeriodic - start_n)* (end_n - start_n)' / l2));
                projection = start_n + t*v1;  % Projection falls on the segment
                distance = sum((ptsPeriodic-projection).^2,2);
            
                
%                 t = max(0, min(1, (pts - start_n)* (end_n - start_n)' / l2));
%                 projection = start_n + t*v1;  % Projection falls on the segment
%                 distance = sum((pts-projection).^2,2);
             end
           % ptsIn(distance<=radius(j)^2) = true;
           
            
        case 3
            % determine alpha and beta are acute angle
            alpha = acosd(((pts - start_n)*(end_n - start_n)')...
                ./(norm(pts - start_n)*norm(end_n - start_n)));
            beta = acosd(((pts - end_n)*(start_n - end_n)')...
                ./(norm(pts - end_n)*norm(start_n - end_n)));
            distance = zeros(n^3,1);
            notAcute = and(alpha<90,beta<90);
            % if not acute angle, distance to line
            distance(notAcute) = vecnorm(cross(repmat(end_n - start_n,sum(notAcute),1),pts(notAcute,:) - start_n),2,2)...
                ./norm(end_n - start_n);
            % if it is acute angle, distance to node
            distance(~notAcute) = min(vecnorm(pts(~notAcute,:) - start_n,2,2),vecnorm(pts(~notAcute,:) - end_n,2,2));
          %  
    end
    %ptsIn = min(distance-radius(j)^2,ptsIn);
    distance = min(reshape(distance,[],9),[],2);
    ptsIn(distance<radius(j)^2+tol) = true;
    thisSignedDist = sqrt(distance)-radius(j);
    signedDist = -max(-thisSignedDist,-signedDist);
end
%ptsIn = 1./(1+exp(ptsIn*10/radForSigmoid));
end