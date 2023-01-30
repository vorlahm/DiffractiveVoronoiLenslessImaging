function isin = isinner(c, lb, ub)
%ISINNER detects if the centers c are inside the boundary [lb, ub].
%   ISIN = ISINNER(C, LB, UB) detects if the centers c are inside of the
%   boundaries defined by LB and UB.
%
%   author: Qiang Fu
%   qiang.fu@kaust.edu.sa
%   2023-01-31

b1 = repmat(lb, [1 size(c, 2)]);
b2 = repmat(ub, [1 size(c, 2)]);
d1 = c - b1;
d2 = c - b2;

if all(d1(:) >= 0) && all(d2(:) <= 0)
    isin = true;
else 
    isin = false;
end
end