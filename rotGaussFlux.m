function [fx, fy] = rotGaussFlux( t, x1, x2, c)
fx = -4.*x2 .* c;
fy =  4.*x1 .* c;
end