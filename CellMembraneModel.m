function CellMembraneModel()
clf;

tInitial = 0;
tFinal = 1;
dt = 0.1;
t = 0:dt:1;

testBodyForce = [1, 0];
restSpringLength = 0.45;
springConst = 5;

maxNodeDist = 1;
minNodeDist = 0.1;

L = linspace(0,2.*pi,9);
L = L(1:end-1);
x = 1.2*cos(L)';
y = 1.2*sin(L)';

r = [x, y];

numNodes = size(r,1);

plotCell(r);
pause(0.5);

for i = 1:length(t)
    
    r_prev = r;
    
    
    
    for currentNode = 1:numNodes
        
        springForce = calcSpringForce(currentNode);
        
        
        r(currentNode,:) = r_prev(currentNode,:) + dt * (springForce);
    end
    r(1,:) = r_prev(1,:) + dt * (testBodyForce);
    
    
    
    %add/remove nodes that get too far apart
    r_new = addNodes(r);
    r = r_new;
    numNodes = size(r,1);
    
    if numNodes > 3
        r_new = removeNodes(r);
        r = r_new;
        numNodes = size(r,1);
    end
    
    plotCell(r);
    title(['frame = ' num2str(i)]);
    pause(0.5);
    
end

    function plotCell(r)
        rPlot = [r(:,:); r(1,:)];
        plot(rPlot(:,1), rPlot(:,2)+.2, 'rx-');
        axis(2.4*[-1 1 -1 1]);
        
        for node = 1:numNodes
            text(r(node,1), r(node,2), num2str(node), 'FontSize', 16);
        end
    end

%%
    function restoreForce = calcRestoreForce(currentnode)
        prevNode = mod(currentNode-1,numNodes) + numNodes*(currentNode == 1);
        nextNode = mod(currentNode+1,numNodes) + numNodes*(currentNode == numNodes-1);
         %%%%
        
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
        
        r_new = r;
        
        for currentNode = 1:size(r,1)
            nextNode = mod(currentNode+1,numNodes) + numNodes*(currentNode == numNodes-1);
            
            dist = sqrt(sum((r(nextNode,:) - r(currentNode,:)).^2));
            
            if dist > maxNodeDist
                
                midpoint = (r(nextNode,:) + r(currentNode,:))/2;
                
                if currentNode == numNodes
                    r_new = [r_new(1:end,:); midpoint];
                else
                    r_new = [r_new(1:currentNode,:); midpoint; r_new(nextNode:end,:)];
            
                end
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
        end
    end
        
        
    
        






end