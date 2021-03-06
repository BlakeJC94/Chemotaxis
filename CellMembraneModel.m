function CellMembraneModel()
% Model of cell membrane
%
% This model takes n ordered points in an anti-clockwise polygon and joins
% neighboring nodes with 'springs' and applys  varior curvature and volume
% conservation forces
%
% Nodes are numbered 0, 1, 2, ..., (numNodes-1) and the connections between
% 0 and 1, 1 and 2, ..., (numNodes-1) and 0 are numbered in the same way.
% Equation of motion overdampened (inertia-free) and solved with Forward
% Euler method with timestep dt
%
% Based on concpets from Phagocytosis model by David M. Richards
% Blake Cook - 26/01/2018
%% Preamble
clf;

%% Set up parameters
tInitial = 0;
tFinal = 4;
dt = 0.05;

testBodyForce = [0.1, 0];


springConst = 2;

maxNodeDist = 1;
minNodeDist = 0.1;

curvConst = 3;

cellVolume = 4;
volConst = 0.5;

%% initialize nodes and connections
L = linspace(0,2.*pi,9);
L = L(1:end-1);
x = 1.2*cos(L)';
y = 1.2*sin(L)';

r = [x, y];
r(1,:) = r(1,:) + [1, 0];

r(4,:) = r(3,:) - [0.3, 0];
r(6,:) = r(7,:) - [0.3, 0];

restSpringLengths = 0.9 * ones(size(r,1),1);
% restSpringLengths = 0.9 * sqrt(sum(([diff(r);r(1,:)-r(end,:)]).^2,2));

plotCell(r);
pause(0.1);

%% run simulation
for t = (tInitial+dt):dt:tFinal
    
    rPrev = r;
    numNodes = size(r,1);
    currentCellVolume = polyarea(r(:,1), r(:,2));
    
    % loop over nodes
    for currentNode = 0:(numNodes-1)
        
        restAngle = pi*(1-2/numNodes);
        
        force = [0, 0];
        
        springForce = calcSpringForce(currentNode);
        force = force + springForce;
        
        curvForce = calcCurvatureForce(currentNode);
        force = force + curvForce;
        
        volForce = calcVolumeForce(currentNode);
        if sum(isnan(volForce))>0
            error('NANS IN volForce');
        end
        force = force + volForce;
        
        force = force + testBodyForce*(currentNode == 0);
        
        if sum(isnan(force))>0
            error('NANS IN force');
        end
        
        r(currentNode+1,:) = rPrev(currentNode+1,:) + dt * force;
        
    end
    
    %add/remove nodes that get too far apart
    [rNew, restSpringLengthsNew] = updateNodes(r);
    r = rNew;
    restSpringLengths = restSpringLengthsNew;
    
    plotCell(r);
    title(['t = ' num2str(t)]);
    pause(0.1);
end

%%
    function plotCell(r)
        rPlot = [r(:,:); r(1,:)];
        plot(rPlot(:,1), rPlot(:,2)+.2, 'rx-');
        axis(2.4*[-1 1 -1 1]);
        
        for node = 0:(size(r,1)-1)
            text(r(node+1,1), r(node+1,2), num2str(node), 'FontSize', 16);
            nextNode = mod(node+1,size(r,1));
            
            midpoint = (r(node+1,:) + r(nextNode+1,:))/2;
            text(midpoint(1), midpoint(2), sprintf('%0.2f',restSpringLengths(node+1)), 'FontSize', 10);
        end
    end

%%
    function volForce = calcVolumeForce(node)
        
        volForce = [0, 0];
        
        prevNode = mod(node-1,numNodes);
        nextNode = mod(node+1,numNodes);
        
        dirPrev = rPrev(prevNode+1,:) - rPrev(node+1,:);
        dirPrev = dirPrev/norm(dirPrev);
        
        dirNext = rPrev(nextNode+1,:) - rPrev(node+1,:);
        dirNext = dirNext/norm(dirNext);
        
        dir = (dirPrev + dirNext);
        if norm(dir) == 0
            dir = [-dirNext(2), dirNext(1)];
        else
            dir = (dirPrev + dirNext)/norm(dirPrev + dirNext);
        end
        
        rExclCurrent = rPrev(~ismember(1:size(rPrev,1), (node+1)),:);
        
        nodeInCell = inpolygon(rPrev(node+1,1),rPrev(node+1,2),...
            rExclCurrent(:,1),rExclCurrent(:,2));
        
        volForce = volForce + volConst * dir * (1-2*nodeInCell) * (currentCellVolume - cellVolume);

    end

%%
    function curveForce = calcCurvatureForce(node)
        
        curveForce = [0, 0];
        
        prevNode = mod(node-1,numNodes);
        nextNode = mod(node+1,numNodes);
        
        currentPoint = rPrev(node+1,:);
        prevPoint = rPrev(prevNode+1,:) - currentPoint;
        nextPoint = rPrev(nextNode+1,:) - currentPoint;
        
        angle = acos(dot(prevPoint, nextPoint)/(norm(prevPoint)*norm(nextPoint)));
        
        % check if currentPoint is in the polygon formed by removing the
        % currentPoint
        rExclCurrent = rPrev(~ismember(1:size(rPrev,1), (node+1)),:);
        
        nodeInCell = inpolygon(currentPoint(1),currentPoint(2),...
            rExclCurrent(:,1),rExclCurrent(:,2));
        
        
        if (angle < restAngle) && (~nodeInCell)
            % half angle between nextPoint and prevPoint
            curveForce = (prevPoint/norm(prevPoint)) + (nextPoint/norm(nextPoint));
            curveForce = curvConst * (restAngle - angle) * curveForce/norm(curveForce);
        end
        
    end

%%
    function springForce = calcSpringForce(node)
        
        springForce = [0, 0];
        
        % calculate spring force from previous node
        prevNode = mod(node-1,numNodes);
        restSpringLength = restSpringLengths(prevNode+1);
        dist = norm(rPrev(prevNode+1,:) - rPrev(node+1,:));
        dir = (rPrev(prevNode+1,:) - rPrev(node+1,:))/dist;
        springForce = springForce + springConst * dir * (dist - restSpringLength);
        
        % calculate spring force from next node
        nextNode = mod(node+1,numNodes);
        restSpringLength = restSpringLengths(node+1);
        dist = norm(rPrev(nextNode+1,:) - rPrev(node+1,:));
        dir = (rPrev(nextNode+1,:) - rPrev(node+1,:))/dist;
        springForce = springForce + springConst * dir * (dist - restSpringLength);
        
        
    end

%%
    function [rNew, restSpringLengthsNew] = updateNodes(r)
        
        % Disallow (3/5)*maxNodeDist <= minNodeDist to prevent adding
        % inifinitely many new nodes (more conservative bound)
        if (3/5)*maxNodeDist <= minNodeDist
            error('(3/5)*maxNodeDist <= minNodeDist, likely will result in infinitely new nodes')
        end
        
        %DEBUG: error if there is NANS in input
        if sum(isnan(r(:)))
            error('NANS IN updateNodes INPUT')
        end
        
        % initialise output
        rNew = r;
        restSpringLengthsNew = restSpringLengths;
        
        
        
        % % NODE REMOVAL
        
        % check lengths and check which nodes need to be removed
        connDists = sqrt(sum(([diff(r);r(1,:)-r(end,:)]).^2,2));
        removeRows = (connDists < minNodeDist);
        
        % if there's any rows to be removed
        if sum(removeRows) > 0
            
            % stop if too many rows are removed
            if numNodes - sum(removeRows) < 4
                error('Too many rows removed, try reducing minNodeDist');
            end
            
            % mark rows for removal
            rNew(removeRows,:) = NaN;
            restSpringLengthsNew(removeRows) = NaN;
            
            % replace removed nodes with midpoints
            if contains(num2str([removeRows', removeRows(1)]), num2str([1 1]))
                % if there are blocks of rows to be removed, identify
                % blocks and average positions over each block
                
                % identify blocks 
                %  - blockstart: first row of ones
                %  - blockend: row after last row of ones in block
                blockstart = strfind(removeRows',[0 1])+1;
                blockend = strfind(removeRows',[1 0])+1;
                % rearrange blockend if there is a block that overlaps end
                if blockend(1) < blockstart(1)
                    blockend = [blockend(2:end), blockend(1)];
                end
                
                for index = 1:length(blockstart)
                    % for each block, get indices of removed rows
                    startrow = blockstart(index);
                    endrow = blockend(index);
                    % (get indices of block that overlaps end of nodes)
                    if (index == length(blockstart)) && (endrow < startrow)
                        endrow = numNodes + endrow;
                    end
                    
                    samplerows = mod((startrow-1):(endrow-1),numNodes)+1;
                    % replace value in row 'blockend' with average over all
                    % the elements in a block
                    rNew(mod(endrow-1,numNodes)+1,:) = sum(r(samplerows,:),1)/numel(samplerows);
                    
                    % replace previous spring of block with average 
                    restSpringLengthsNew(mod(removeRows-2,numNodes)+1) = ...
                        sum(restSpringLengths(samplerows,:),1)/numel(samplerows);
                    
                end
                
            else
                % if only isolated rows need to be removed, relpace the
                % value in the next row with the average
                removeRows = find(removeRows);
                rNew(mod(removeRows,numNodes)+1,:) = ...
                    (r(mod(removeRows,numNodes)+1,:)+r(removeRows,:))/2;
                
                % replace previous springs of removed nodes
                restSpringLengthsNew(mod(removeRows-2,numNodes)+1) = ...
                    (restSpringLengths(mod(removeRows-2,numNodes)+1) ...
                    + restSpringLengths(mod(removeRows-1,numNodes)+1))/2;
                
            end
            
            % remove marked nodes
            rNew(isnan(rNew(:,1)),:) = [];
            restSpringLengthsNew(isnan(restSpringLengthsNew(:,1))) = [];
            % update connection distances with removed/replaced nodes
            connDists = sqrt(sum(([diff(rNew);rNew(1,:)-rNew(end,:)]).^2,2));
            
        end
        
        
        
        % % ADD NODES
        
        % check which nodes need to be added
        addRows = (connDists > maxNodeDist);
        
        if sum(addRows) > 0
            % if a row need to be added
            
            % update existing spring lengths
            restSpringLengthsNew(addRows) = restSpringLengthsNew(addRows)...
                - connDists(addRows)/2;
            
            % preallocate new positions and new spring lengths
            numNewRows = sum(addRows);
            rNew_tmp = NaN*ones(size(rNew,1)+numNewRows,2);
            restSpringLengthsNew_tmp = NaN*ones(size(rNew,1)+numNewRows,1);
            
            % carry over original points, leave NaN rows for new nodes
            indices = (1:size(rNew,1)) + cumsum([0, addRows(1:end-1)']);
            rNew_tmp(indices,:) = rNew(:,:);
            restSpringLengthsNew_tmp(indices) = restSpringLengthsNew(:);
            
            
            
            % fill neach NaN row with the midpoint and prev springlength
%             addRows = find(addRows);
%             for index = 1:numNewRows
%                 row = addRows(index);
%                 nextrow = mod(row,size(rNew,1))+1;
%             
%                 rNew_tmp(row+index,:) = (rNew(row,:)+rNew(nextrow,:))/2;
%                 restSpringLengthsNew_tmp(row+index) = restSpringLengthsNew_tmp(row+index-1);
%             end
            
            addRows = find(addRows);
            indices = (1:length(addRows))' + addRows;
            addRowsNext = mod(addRows,size(rNew,1))+1;
            rNew_tmp(indices,:) = (rNew(addRows,:)+rNew(addRowsNext,:))/2;
            restSpringLengthsNew_tmp(indices) = restSpringLengthsNew_tmp(indices-1);

            rNew = rNew_tmp;
            restSpringLengthsNew = restSpringLengthsNew_tmp;
            
        end
        
        
        %DEBUG: error if there is NANS in output
        if sum(isnan(rNew(:)))
            error('NANS IN updateNodes INPUT')
        end
        
        
    end


%%
%     function [rNew, restSpringLengthsNew] = updateNodes(r)
%
%         % Disallow (3/5)*maxNodeDist <= minNodeDist to prevent adding
%         % inifinitely many new nodes (more conservative bound)
%         if (3/5)*maxNodeDist <= minNodeDist
%             error('(3/5)*maxNodeDist <= minNodeDist, likely will result in infinitely new nodes')
%         end
%
%         rNew = r;
%         restSpringLengthsNew = restSpringLengths;
%         %%%%
%
%         % col 1: -1 if removed, +1 if new point added after
%         pointsUpdate = zeros(numNodes, 3);
%
%
%         for node = 0:(numNodes-1)
%
%             nextNode = mod(node+1, numNodes);
%             prevNode = mod(node-1, numNodes);
%             dist = sqrt(sum((rNew(nextNode+1,:) - rNew(node+1,:)).^2));
%
%             if dist < minNodeDist
%                 % if connection is too small, remove node and replace
%                 % nextNode with midpoint between node and nextNode (so
%                 % search direction doesnt influence results)
%                 % can be visualized as node and nextNode merging
%                 pointsUpdate(node+1,1) = -1;
%
%                 midpoint = (rNew(nextNode+1,:) + rNew(node+1,:))/2;
%                 rNew(nextNode+1,:) = midpoint;
%
%                 rNew(node+1,:) = [NaN, NaN];
%                 restSpringLengthsNew(node+1) = NaN;
%
%
%             elseif dist > maxNodeDist
%                 % if connection is too large, add new point after node at
%                 % midpoint between node and nextNode
%                 pointsUpdate(node+1,1) = +1;
%
%                 restSpringLengthsNew(node+1) = restSpringLengthsNew(node+1) - dist/2;
%
%                 midpoint = (rNew(nextNode+1,:) + rNew(node+1,:))/2;
%                 pointsUpdate(node+1,2:3) = midpoint;
%
%
%
%             end
%
%         end
%
%         % check if too many points are being removed
%         if length(find(pointsUpdate(:,1)>-1)) < 4
%             error('Too many nodes removed, try decreasing minNodeDist');
%         end
%
%         % add nodes
%         [rowsUpdate, ~] = find(pointsUpdate(:,1)==+1);
%
%         pointsUpdate = [rowsUpdate, pointsUpdate(rowsUpdate,2:3)];
%
%         numNewRows = numel(rowsUpdate);
%         rNew = [rNew; NaN*ones(numNewRows,2)];
%         restSpringLengthsNew = [restSpringLengthsNew; NaN*ones(numNewRows,1)];
%         for index = 1:length(rowsUpdate)
%
%             row = pointsUpdate(index,1)+index-1;
%
%             rNew = [rNew(1:row,:);...
%                 pointsUpdate(index, 2:3);...
%                 rNew((row+1):(end-numNewRows+index),:)];
%
%             restSpringLengthsNew = [restSpringLengthsNew(1:row);...
%                 restSpringLengthsNew(row);...
%                 restSpringLengthsNew((row+1):(end-numNewRows+index))];
%
%         end
%
%         % remove nodes
%         rNew(isnan(rNew(:,1)),:) = [];
%         restSpringLengthsNew(isnan(restSpringLengthsNew(:,1))) = [];
%     end




end