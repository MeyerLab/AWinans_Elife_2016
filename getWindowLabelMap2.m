function [labelMask, boundaryOut]=getWindowLabelMap2(cellMask,windowCoors)
% This function creates a label matrix mapping points near the edge of the
% cell to the closest parametrized edge point (specified in windowCoors).

%written by S. Collins

% ringMask=cellMask-imerode(cellMask,strel('disk',edgeDepthDist));
edgeCoorMask=false(size(cellMask));
edgeCoorPixInd=windowCoors(:,1)+(windowCoors(:,2)-1)*size(cellMask,1);
edgeCoorMask(edgeCoorPixInd)=true;
[D,L]=bwdist(edgeCoorMask);
labelMask=L;
for i=1:size(windowCoors,1)
    labelMask(L==edgeCoorPixInd(i))=i;
end

perimeterAll = false(size(cellMask));
for i = 1:size(windowCoors, 1)
    
    subMask = labelMask == i;
    perimeter = bwperim(subMask);
    perimeterAll = perimeterAll + perimeter;
end
    
[x, y] = find(perimeterAll);
edgeInd = x+(y-1)*size(cellMask,1);

for i = 1:length(edgeInd)
    edgeInd(i) = cellMask(edgeInd(i));
end


x(edgeInd == 0) = [];
y(edgeInd == 0) = [];
boundaryOut = [x y];

labelMask(cellMask==0)=0;
