function isb = isboundary(v, lb, ub)
%ISBOUNDARY detects if the vertices v are on the boundary.
%   ISB = ISBOUNDARY(V, LB, UB) detects if the vertices v are on the
%   boundaries defined by LB and UB.
%
%   author: Qiang Fu
%   qiang.fu@kaust.edu.sa
%   2023-01-31

b1 = repmat(lb, [1 size(v, 2)]);
b2 = repmat(ub, [1 size(v, 2)]);
d1 = v - b1;
d2 = v - b2;

if any(d1(:) == 0) || any(d2(:) == 0)
    isb = true;
else 
    isb = false;
end
end