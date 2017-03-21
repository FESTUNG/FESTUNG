function ret = generateFlipArray(g)
K = g.numT;
ret = true( K, 3 );
for n=1:3
    ret(:, n) = g.T0E( g.E0T(:, n), 2) ~= (1:K)';
end

end % function