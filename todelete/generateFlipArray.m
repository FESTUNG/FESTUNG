function ret = generateMarkSideE0T(g)
K = g.numT;
ret = true( K, 3, 2 );
for n=1:3
    ret(:, n, 1) = g.T0E( g.E0T(:, n), 2) ~= (1:K)';
    ret(:, n, 2) = ~ret(:, n, 2);
end
end % function