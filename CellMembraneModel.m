function CellMembraneModel()

tInitial = 0;
tFinal = 1;
dt = 0.1;
t = 0:dt:1;

testBodyForce = [1, 0];
restSpringLength = 0.8;
springConst = 20;

L = linspace(0,2.*pi,9);
L = L(1:end-1);
x = 1.2*cos(L)';
y = 1.2*sin(L)';

r = [x, y];

numNodes = size(r,1);

rPlot = [r(:,:); r(1,:)];
plot(rPlot(:,1), rPlot(:,2)); 
axis(1.2*[-1 1 -1 1]);
pause(0.5);

for i = 1:length(t)
    r_prev = r;
    
    for currentNode = 1:numNodes
        
        
        nextNode = mod(currentNode+1,numNodes) + numNodes*(currentNode == numNodes-1);
        prevNode = mod(currentNode-1,numNodes) + numNodes*(currentNode == 1);
        
        
        distance = sqrt(sum((r_prev(nextNode,:) - r_prev(currentNode,:)).^2));
        dir = (r_prev(nextNode,:) - r_prev(currentNode,:))/distance;
        springForce = dir * springConst * (distance - restSpringLength);
        
        distance = sqrt(sum((r_prev(prevNode,:) - r_prev(currentNode,:)).^2));
        dir = (r_prev(prevNode,:) - r_prev(currentNode,:))/distance;
        springForce = springForce + dir * springConst * (distance - restSpringLength);
        
        
        bodyForce = 0;
        r(currentNode,:) = r_prev(currentNode,:) + dt * (springForce);
    end
    
    
%     r(1,:) = r_prev(1,:) + dt * (testBodyForce);
    
    rPlot = [r(:,:); r(1,:)];
    plot(rPlot(:,1), rPlot(:,2)); 
    axis(1.2*[-1 1 -1 1]);
    title(['frame = ' num2str(i)]);
    pause(0.5);
    
end



    function calcSpringForce



end