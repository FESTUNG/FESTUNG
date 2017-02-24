function [fx, fy] = evalLinearAdvectionFlux( t, x1, x2, c)
assert(isequal(size(x1), size(x2)), 'x1 and x2 must have the same size')
fx = 1 .* c;
fy = 1 .* c;
end