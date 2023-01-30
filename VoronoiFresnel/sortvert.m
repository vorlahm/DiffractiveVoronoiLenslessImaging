function sortedvert = sortvert(vert)
% SORTVERT sorts the unordered vertices to form a proper polygon.
%
%   SORTEDVERT = SORTVERT(VERT) sorts the input unordered vertices in order
%   to form a proper convex polygon. The center is defined as the mean of 
%   x and y coordinates of all the vertices. The vertices are then sorted 
%   by the angles formed by the line that connects the vertex and the 
%   center with respect to the horizon. The sorted vertices are returned as
%   sortedvert.
%   
%   vert: An unsorted list of vertices.
%
%   author: Qiang Fu
%   qiang.fu@kaust.edu.sa
%   2023-01-31

x = vert(1, :);
y = vert(2, :);
meanx = mean(x);
meany = mean(y);
angles = atan2(y - meany, x - meanx);
[~, idx] = sort(angles);
reorderedx = x(idx);
reorderedy = y(idx);
sortedvert = [reorderedx; reorderedy];
end