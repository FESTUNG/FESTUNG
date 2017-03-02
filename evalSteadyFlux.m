function [fx, fy] = evalSteadyFlux( t, x1, x2, c)
assert(isequal(size(x1), size(x2)), 'x1 and x2 must have the same size')
fx = exp((x1+x2)/2) .* c;
fy = exp((x1-x2)/2) .* c;
end