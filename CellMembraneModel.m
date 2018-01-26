function CellMembraneModel()
clf;

tInitial = 0;
tFinal = 1;
dt = 0.1;
t = 0:dt:1;

testBodyForce = [1, 0];
restSpringLength = 0.9;
springConst = 2;

maxNodeDist = 1;
minNodeDist = 0.1;

% restAngle = 3*pi/4; %angle between points in octagon (set to pi*(1-2/length(r)))
restoreConst = 1;

L = linspace(0,2.*pi,9);
L = L(1:end-1);
x = 1.2*cos(L)';
y = 1.2*sin(L)';

r = [x, y];

r(1,:) = r(1,:) + [1, 0];

numNodes = size(r,1);

plotCell(r);
pause(0.5);

for i = 1:length(t)
    
    r_prev = r;

    for currentNode = 1:numNodes
        
        restAngle = pi*(1-2/numNodes);
        
        springForce = calcSpringForce(currentNode);
        restoreForce = calcRestoreForce(currentNode);

        force = springForce + restoreForce;
        r(currentNode,:) = r_prev(currentNode,:) + dt * (force);
    end
%     r(1,:) = r_prev(1,:) + dt * (testBodyForce);
    
%     testBodyForce = calcRestoreForce(1);
%     r(1,:) = r_prev(1,:) + dt * (testBodyForce);
    
    
    %add/remove nodes that get too far apart
    r_new = addNodes(r);
    r = r_new;
    numNodes = size(r,1);
    
    if numNodes > 4
        r_new = removeNodes(r);
        r = r_new;
        numNodes = size(r,1);
    end
    
    plotCell(r);
    title(['frame = ' num2str(i)]);
    pause(0.5);
    
end

%%
    function plotCell(r)
        rPlot = [r(:,:); r(1,:)];
        plot(rPlot(:,1), rPlot(:,2)+.2, 'rx-');
        axis(2.4*[-1 1 -1 1]);
        
        for node = 1:numNodes
            text(r(node,1), r(node,2), num2str(node), 'FontSize', 16);
        end
    end

%%
    function restoreForce = calcRestoreForce(currentNode)
        restoreForce = 0;
        
        prevNode = mod(currentNode-1,numNodes) + numNodes*(currentNode == 1);
        nextNode = mod(currentNode+1,numNodes) + numNodes*(currentNode == numNodes-1);
        
        prevPoint = r(prevNode,:) - r(currentNode,:);
        nextPoint = r(nextNode,:) - r(currentNode,:);
        
        angle = acos(dot(prevPoint, nextPoint)/(norm(prevPoint)*norm(nextPoint)));
        
        if angle < restAngle
            bridge = nextPoint - prevPoint;
            restoreForce = prevPoint - dot(prevPoint,bridge) * bridge/norm(bridge);
            restoreForce = (restoreForce/norm(restoreForce)) * restoreConst;
        end
        
    end

%%
    function springForce = calcSpringForce(currentNode)
        prevNode = mod(currentNode-1,numNodes) + numNodes*(currentNode == 1);
        nextNode = mod(currentNode+1,numNodes) + numNodes*(currentNode == numNodes-1);
        
        dists = sqrt([sum((r_prev(prevNode,:) - r_prev(currentNode,:)).^2), ...
            sum((r_prev(nextNode,:) - r_prev(currentNode,:)).^2)]);
        
        dir = [(r_prev(prevNode,:) - r_prev(currentNode,:))/dists(1);...
            (r_prev(nextNode,:) - r_prev(currentNode,:))/dists(2)];
        
        springForce = springConst * ( dir(1,:) * (dists(1) - restSpringLength) + ...
            dir(2,:) * (dists(2) - restSpringLength) );
        
    end

%%
    function r_new = addNodes(r)
        % compare currentNode and nextNode and add a new node in midpoint if
        % they're too far apart
        
        newPoints = [];
        
        for currentNode = 1:size(r,1)
            nextNode = mod(currentNode+1,numNodes) + numNodes*(currentNode == numNodes-1);
            dist = norm(r(nextNode,:) - r(currentNode,:));
            
            if dist > maxNodeDist
                
                midpoint = (r(nextNode,:) + r(currentNode,:))/2;
                newPoints = [newPoints; currentNode, midpoint];
                
            end
        end
        
        r_new = [r; zeros(size(newPoints,1),2)];
        for newRow = 1:size(newPoints,1)
            newNode = newPoints(newRow,1);
            midpoint = newPoints(newRow,2:3);
            if newNode == numNodes
                r_new = [r_new(1:end,:); midpoint];
            else
                r_new = [r_new(1:newNode,:); midpoint; r(newNode+1:end,:)];
            end
        end
    end

%%
    function r_new = removeNodes(r) 
        % compare curentNode and nextNode as remove currentNode if too
        % close together
        
        r_new = r;
        
        for currentNode = 1:size(r,1)
            nextNode = mod(currentNode+1,numNodes) + numNodes*(currentNode == numNodes-1);
            
            dist = sqrt(sum((r(nextNode,:) - r(currentNode,:)).^2));
            
            if dist < minNodeDist
                disp(currentNode);
                disp('TRACE');
                
                if currentNode == 1
                    r_new = [r_new(2:end,:)];
                elseif currentNode == numNodes
                    r_new = [r_new(1:end-1,:)];
                else
                    r_new = [r_new(1:currentNode-1,:); r_new(currentNode+1:end,:)];
                end
            
            end
            
            if size(r_new,1) < 4
                error('Too many nodes removed, try decreasing minNodeDist');
            end
            
        end
        
    end
        
        
    
        






end