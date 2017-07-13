function [fx, fy] = evalLeVequeFlux( t, x1, x2, c)
assert(isequal(size(x1), size(x2)), 'x1 and x2 must have the same size')
fx = (0.5 - x2) .* c;
fy = (x1 - 0.5) .* c;
end